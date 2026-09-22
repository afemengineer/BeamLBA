import json

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

from beamlba.assessment import assess
from beamlba.benchmarks import crowded_frame, fixture_metadata, vlasov_beam
from beamlba.bundle import load_bundle, pack_matrix_market


def test_sparse_native_vectors_matrix_market(tmp_path):
    pencil, _ = vlasov_beam(4)
    values, vectors = pencil.dense()
    mmwrite(tmp_path/'k.mtx', pencil.k)
    mmwrite(tmp_path/'g.mtx', pencil.g)
    mmwrite(tmp_path/'native.mtx', sparse.coo_matrix(vectors[:, :3]))
    spec = {'metadata': fixture_metadata('sparse native vectors'),
            'members': [{'name': 'M', 'dofs': list(range(pencil.n))}],
            'native_vectors_mtx': 'native.mtx', 'native_values': values[:3].tolist()}
    (tmp_path/'spec.json').write_text(json.dumps(spec))
    pack_matrix_market(str(tmp_path/'k.mtx'), str(tmp_path/'g.mtx'),
                       str(tmp_path/'spec.json'), str(tmp_path/'bundle'))
    loaded, _, native = load_bundle(tmp_path/'bundle')
    assert loaded.native_parity(*native)['maximum_residual'] < 1e-8


def test_clustered_local_seeds_preserve_targeting():
    pencil = crowded_frame(30)
    pencil.members = pencil.members[20:25]
    pencil.config.max_shifts = 5
    report, vectors = assess(pencil)
    assert all(row['candidate_modes'] for row in report['members'])
    assert all(pencil.residual(row['factor'], vectors[:, i]) < 1e-7
               for i, row in enumerate(report['modes']))
    assert len(report['searches']) <= 5
