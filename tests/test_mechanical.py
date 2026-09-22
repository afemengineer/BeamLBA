import importlib.util
from pathlib import Path
from types import SimpleNamespace as NS

import pytest

spec = importlib.util.spec_from_file_location('mechanical_launcher', Path(__file__).parents[1]/'mechanical'/'beamlba_mechanical.py')
launcher = importlib.util.module_from_spec(spec)
spec.loader.exec_module(launcher)


def test_mesh_element_selection():
    selection = NS(Name='MEMBER_A', Location=NS(SelectionType='MeshElements', Ids=[3, 1, 3]))
    assert launcher.member_elements(selection, None, None) == [1, 3]


def test_line_body_selection():
    selection = NS(Name='MEMBER_A', Location=NS(SelectionType='GeometryEntities', Ids=[12]))
    mesh = NS(MeshRegionById=lambda _: NS(ElementIds=[1, 2]))
    geo = NS(GeoEntityById=lambda _: NS(BodyType='GeoBodyWire'))
    assert launcher.member_elements(selection, mesh, geo) == [1, 2]


@pytest.mark.parametrize('kind', ['MeshNodes', 'MeshElementFaces'])
def test_nodal_or_face_ownership_rejected(kind):
    with pytest.raises(ValueError):
        launcher.member_elements(NS(Name='M', Location=NS(SelectionType=kind, Ids=[1])), None, None)


def test_surface_is_not_line_body():
    selection = NS(Name='M', Location=NS(SelectionType='GeometryEntities', Ids=[12]))
    geo = NS(GeoEntityById=lambda _: NS(BodyType='GeoBodySheet'))
    with pytest.raises(ValueError):
        launcher.member_elements(selection, None, geo)


def test_analysis_ambiguity():
    with pytest.raises(ValueError):
        launcher.select_analysis([NS(Name='LBA'), NS(Name='LBA')], 'LBA')


def test_worker_arguments_no_shell(tmp_path):
    args = launcher.worker_arguments(str(tmp_path/'model with spaces'), str(tmp_path/'out'))
    assert args[:3] == ['-m', 'beamlba.cli', 'run']
    assert len(args) == 6
    assert launcher.RUN_ASSESSMENT is False


def test_collect_read_only_snapshot():
    selected = NS(Name='MEMBER_A', ObjectId=101, Location=NS(SelectionType='MeshElements', Ids=[1]), Children=[])
    model = NS(Analyses=[NS(Name='Static', ObjectId=1, AnalysisType='StaticStructural'),
                         NS(Name='LBA', ObjectId=2, AnalysisType='EigenvalueBuckling', Solution=NS(Children=[]))],
               NamedSelections=NS(Children=[selected]))
    mesh = NS(ElementById=lambda _: NS(NodeIds=[1, 2], Type='Line2'),
              NodeById=lambda i: NS(X=float(i), Y=0., Z=0.))
    data = NS(Project=NS(Model=model), MeshDataByName=lambda _: mesh, GeoData=None)
    inventory = launcher.collect(NS(DataModel=data), 'Static', 'LBA', 'MEMBER_')
    assert inventory['status'] == 'INVENTORY_ONLY'
    assert inventory['snapshot']['members'][0]['element_ids'] == [1]
    assert inventory['solver_equation_mapping_available'] is False
    assert len(inventory['snapshot_sha256']) == 64
