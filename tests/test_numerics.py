import numpy as np
import pytest
from scipy import linalg, sparse

from beamlba.assessment import assess
from beamlba.benchmarks import crowded_frame, uniform_mcr, vlasov_beam
from beamlba.numerics import Config, InvalidModel, Member, Pencil, positive_roots


def test_indefinite_and_singular_g():
    values, vectors = positive_roots(np.eye(3), np.diag([-0.5, 0.2, 0]))
    assert values == pytest.approx([2.0])
    assert np.linalg.norm((np.eye(3) + values[0]*np.diag([-0.5, 0.2, 0])) @ vectors[:, 0]) < 1e-12


@pytest.mark.parametrize('elements', [2, 4, 8, 16])
def test_vlasov_convergence(elements):
    pencil, analytical = vlasov_beam(elements)
    values, vectors = pencil.dense()
    error = abs(values[0]*20000/analytical - 1)
    assert error < 0.01
    assert pencil.residual(values[0], vectors[:, 0]) < 1e-8
    if elements > 2:
        coarse, _ = vlasov_beam(elements//2)
        coarse_values, _ = coarse.dense()
        assert error < abs(coarse_values[0]*20000/analytical - 1)


def test_reference_16_element_accuracy():
    pencil, analytical = vlasov_beam()
    values, _ = pencil.dense()
    assert abs(values[0]*20000/analytical - 1) < 2e-6
    assert analytical == pytest.approx(112887.32064731122)


def test_predictor_matches_single_beam():
    pencil, _ = vlasov_beam(8)
    prediction, seed = pencil.predict(pencil.members[0])
    values, _ = pencil.dense()
    assert prediction == pytest.approx(values[0], rel=1e-8)
    near, vectors, _ = pencil.near(prediction*(1-1e-5), seed)
    assert min(abs(near-values[0])) < 1e-6
    assert max(pencil.residual(v, x) for v, x in zip(near, vectors.T)) < 1e-6


def test_rank_281_without_preceding_modes():
    pencil = crowded_frame()
    all_values, all_vectors = pencil.dense()
    target = pencil.members[280]
    assert np.max(np.sum(all_vectors[target.dofs, :20]**2, axis=0)) < 1e-8
    pencil.members = [target]
    report, vectors = assess(pencil)
    assert len(report['searches']) == 1
    assert len(report['modes']) <= pencil.config.modes_per_shift
    factors = np.array([m['factor'] for m in report['modes']])
    j = int(np.argmin(abs(factors-all_values[280])))
    assert factors[j] == pytest.approx(all_values[280], rel=1e-10)
    assert pencil.residual(factors[j], vectors[:, j]) < 1e-10
    assert report['members'][0]['status'] == 'CANDIDATE_REQUIRES_REVIEW'
    assert report['coverage_certified'] is False
    assert report['members'][0]['Mcr'] is None


def test_scaled_coordinates_preserve_eigenvalues():
    pencil, _ = vlasov_beam(4)
    values, _ = pencil.dense()
    factors = np.geomspace(1e-4, 1e4, pencil.n)
    d = sparse.diags(factors)
    transformed = Pencil(d @ pencil.k @ d, d @ pencil.g @ d, [])
    transformed_values, _ = transformed.dense()
    assert transformed_values[:3] == pytest.approx(values[:3], rel=1e-8)


def test_native_modes_reused_without_predictors():
    pencil, _ = vlasov_beam(4)
    values, vectors = pencil.dense()
    report, _ = assess(pencil, initial_modes=(values[:3], vectors[:, :3]))
    assert report['initial_mode_count'] == 3
    assert not report['searches']
    assert pencil._elastic_factor is None


def test_native_sign_or_mapping_mismatch_rejected():
    pencil, _ = vlasov_beam(4)
    values, vectors = pencil.dense()
    with pytest.raises(InvalidModel, match='parity failed'):
        pencil.native_parity(-values[:2], vectors[:, :2])
    with pytest.raises(InvalidModel, match='parity failed'):
        pencil.native_parity(values[:2], vectors[::-1, :2])


def test_no_geometric_data_is_unresolved():
    pencil = Pencil(np.eye(4), -np.eye(4), [Member('M', np.array([0, 1]))])
    report, vectors = assess(pencil)
    assert report['members'][0]['status'] == 'UNRESOLVED'
    assert vectors.shape == (4, 0)
    assert report['design_check'] == 'NOT_PERFORMED'


def test_budget_is_obeyed():
    pencil = crowded_frame(20)
    pencil.config.max_shifts = 1
    report, _ = assess(pencil)
    assert len(report['searches']) == 1
    assert any(row['status'] == 'UNRESOLVED' for row in report['members'])


def test_member_size_limit():
    pencil, _ = vlasov_beam(4)
    pencil.config.max_member_dofs = 2
    with pytest.raises(InvalidModel, match='budget'):
        pencil.predict(pencil.members[0])


def test_factorization_is_reused():
    pencil = crowded_frame(6)
    pencil.predict(pencil.members[0])
    factor = pencil._elastic_factor
    pencil.predict(pencil.members[1])
    assert pencil._elastic_factor is factor


@pytest.mark.parametrize('shift', [0, -1, np.nan, np.inf])
def test_invalid_shift(shift):
    with pytest.raises(InvalidModel):
        Pencil(np.eye(4), -np.eye(4), []).near(shift)


@pytest.mark.parametrize('k,g', [
    (np.array([[1, 1], [0, 1]]), -np.eye(2)),
    (np.diag([1, -1]), -np.eye(2)),
    (np.array([[1, 2], [2, 1]]), -np.eye(2)),
    (np.eye(2, dtype=complex)*(1+1j), -np.eye(2)),
    (np.diag([1, np.nan]), -np.eye(2)),
    (np.eye(2), -np.eye(3)),
    (np.zeros((2, 2)), -np.eye(2)),
])
def test_invalid_matrices(k, g):
    with pytest.raises(InvalidModel):
        Pencil(k, g, [])


@pytest.mark.parametrize('dofs', [np.array([0, 0]), np.array([-1]), np.array([4]), np.array([0.5]), np.array([])])
def test_invalid_member_dofs(dofs):
    with pytest.raises(InvalidModel):
        Pencil(np.eye(4), -np.eye(4), [Member('M', dofs)])


def test_overlapping_elements_rejected_shared_nodes_allowed():
    members = [Member('A', np.array([0, 1]), elements=[1]), Member('B', np.array([1, 2]), elements=[1])]
    with pytest.raises(InvalidModel, match='ownership'):
        Pencil(np.eye(4), -np.eye(4), members)
    members[1].elements = [2]
    Pencil(np.eye(4), -np.eye(4), members)


def test_negative_member_energy_rejected():
    with pytest.raises(InvalidModel, match='semidefinite'):
        Pencil(np.eye(4), -np.eye(4), [Member('M', np.array([0, 1]), ke=-np.eye(2))])


@pytest.mark.parametrize('kwargs', [{'max_shifts': 0}, {'max_seconds': -1}, {'cluster_ratio': 1},
                                    {'modes_per_shift': True}, {'participation_min': 2}])
def test_invalid_config(kwargs):
    with pytest.raises(InvalidModel):
        Config(**kwargs).validate()


@pytest.mark.parametrize('iw,length', [(-1, 6), (1, 0), (np.nan, 6)])
def test_invalid_analytical_input(iw, length):
    with pytest.raises(ValueError):
        uniform_mcr(210e9, 8e-6, 80e9, 2e-7, iw, length)


def test_repeated_independent_native_vectors_retained():
    pencil = Pencil(np.eye(4), -np.eye(4), [Member('M', np.array([0, 1]), ke=np.eye(2))])
    report, vectors = assess(pencil, initial_modes=(np.ones(2), np.eye(4)[:, :2]))
    assert len(report['modes']) == 2
    assert np.linalg.matrix_rank(vectors) == 2


def test_unknown_override_rejected():
    with pytest.raises(InvalidModel, match='unknown'):
        assess(Pencil(np.eye(4), -np.eye(4), []), {'absent': 2.0})
