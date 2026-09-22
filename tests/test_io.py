import json

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from beamlba.assessment import assess
from beamlba.benchmarks import fixture_metadata, vlasov_beam
from beamlba.bundle import load_bundle, pack_matrix_market, read_matrix_market, safe_path, write_bundle
from beamlba.cli import main
from beamlba.numerics import InvalidModel
from beamlba.report import write_report


def fixture(tmp_path, source='synthetic', native=False):
    pencil, _ = vlasov_beam(4)
    metadata = fixture_metadata('fixture')
    metadata['source'] = source
    values, vectors = pencil.dense()
    data = (values[:3], vectors[:, :3]) if native else None
    write_bundle(tmp_path/'bundle', pencil.k, pencil.g, pencil.members, metadata, data)
    return tmp_path/'bundle', pencil


def test_roundtrip(tmp_path):
    root, before = fixture(tmp_path, native=True)
    after, manifest, native = load_bundle(root)
    assert np.allclose(before.k.toarray(), after.k.toarray())
    assert after.native_parity(*native)['algebraic_parity']
    assert manifest['metadata']['source'] == 'synthetic'


def test_tamper_rejected(tmp_path):
    root, _ = fixture(tmp_path)
    with (root/'G.npz').open('ab') as stream:
        stream.write(b'tamper')
    with pytest.raises(InvalidModel, match='Integrity'):
        load_bundle(root)


def test_path_escape_rejected(tmp_path):
    with pytest.raises(InvalidModel, match='escapes'):
        safe_path(tmp_path, '../outside.npz')


def test_output_is_not_overwritten(tmp_path):
    root, pencil = fixture(tmp_path)
    with pytest.raises(FileExistsError):
        write_bundle(root, pencil.k, pencil.g, pencil.members, fixture_metadata('fixture'))


def test_triangular_import_is_explicit(tmp_path):
    a = np.array([[2., 1.], [1., 3.]])
    mmwrite(tmp_path/'matrix.mtx', sparse.coo_matrix(np.triu(a)), symmetry='general')
    assert np.allclose(read_matrix_market(tmp_path/'matrix.mtx', 'upper').toarray(), a)
    mmwrite(tmp_path/'full.mtx', a, symmetry='general')
    with pytest.raises(InvalidModel, match='double'):
        read_matrix_market(tmp_path/'full.mtx', 'upper')


def test_ansys_gate_and_diagnostic_status(tmp_path):
    root, _ = fixture(tmp_path, source='ansys')
    assert main(['run', str(root), '--out', str(tmp_path/'blocked')]) == 2
    assert not (tmp_path/'blocked').exists()
    assert main(['run', str(root), '--out', str(tmp_path/'debug'), '--diagnostic']) == 0
    report = json.loads((tmp_path/'debug'/'report.json').read_text())
    assert report['members'][0]['status'] == 'DIAGNOSTIC_ONLY'
    assert not report['input_review_gate_satisfied']


def test_reviewed_ansys_input_still_not_design_pass(tmp_path):
    root, _ = fixture(tmp_path, source='ansys', native=True)
    path = root/'manifest.json'
    manifest = json.loads(path.read_text())
    manifest['metadata'].update(mapping_reviewed=True, element_formulation_reviewed=True)
    path.write_text(json.dumps(manifest))
    assert main(['run', str(root), '--out', str(tmp_path/'report')]) == 0
    report = json.loads((tmp_path/'report'/'report.json').read_text())
    assert report['input_review_gate_satisfied']
    assert report['design_check'] == 'NOT_PERFORMED'
    assert not report['coverage_certified']


def test_report_escapes_names(tmp_path):
    pencil, _ = vlasov_beam(4)
    pencil.members[0].name = '<script>alert(1)</script>'
    report, vectors = assess(pencil)
    write_report(tmp_path/'report', report, vectors)
    html = (tmp_path/'report'/'report.html').read_text()
    assert '<script>' not in html
    assert '&lt;script&gt;' in html


def test_demo_and_verify(tmp_path):
    assert main(['demo', '--kind', 'beam', '--out', str(tmp_path/'demo')]) == 0
    assert main(['verify', str(tmp_path/'demo'/'bundle')]) == 0
    assert (tmp_path/'demo'/'assessment'/'report.html').exists()


def test_nonproportional_input_rejected(tmp_path):
    root, _ = fixture(tmp_path)
    path = root/'manifest.json'
    data = json.loads(path.read_text())
    data['metadata']['load_scaling'] = 'fixed gravity plus live load'
    path.write_text(json.dumps(data))
    with pytest.raises(InvalidModel, match='proportional'):
        load_bundle(root)
