import numpy as np
import pytest
from scipy import sparse

from beamlba.assembly import ElementContribution, assemble_members, assembly_residuals
from beamlba.numerics import InvalidModel


def test_nontrivial_constraint_transform_and_shared_node():
    a = ElementContribution(1, np.array([0, 1]), np.array([[1., 0], [0.5, 0.5], [0, 1.]]), np.eye(3), -np.eye(3))
    b = ElementContribution(2, np.array([1, 2]), np.eye(2), np.eye(2), -2*np.eye(2))
    members = assemble_members({'A': [1], 'B': [2]}, [a, b])
    assert members[0].ke == pytest.approx(a.transform.T @ a.transform)
    k, g = np.zeros((3, 3)), np.zeros((3, 3))
    for member in members:
        k[np.ix_(member.dofs, member.dofs)] += member.ke
        g[np.ix_(member.dofs, member.dofs)] += member.kg
    residuals = assembly_residuals(sparse.csc_matrix(k), sparse.csc_matrix(g), [a, b])
    assert max(residuals.values()) < 1e-14
    bad = assembly_residuals(sparse.csc_matrix(k), sparse.csc_matrix(g), [a])
    assert bad['K_relative_residual'] > 0.1


@pytest.mark.parametrize('selection', [{'A': [1], 'B': [1]}, {'A': [2]}, {'A': []}])
def test_bad_ownership(selection):
    element = ElementContribution(1, np.array([0, 1]), np.eye(2), np.eye(2), -np.eye(2))
    with pytest.raises(InvalidModel):
        assemble_members(selection, [element])


def test_invalid_transform():
    element = ElementContribution(1, np.array([0, 1]), np.eye(3), np.eye(2), -np.eye(2))
    with pytest.raises(InvalidModel):
        element.reduced()
