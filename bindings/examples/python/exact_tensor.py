"""Build an exact rational matrix from row-major cells."""

from fractions import Fraction

import rssn


with rssn.Session() as session:
    wide = (1 << 64) + 1
    with rssn.Value.tensor_cells(session, [wide, Fraction(2, 3), 3, 4], (2, 2)) as matrix:
        with session.call_value("num_native_qmatrix_new", matrix, 1) as normalized:
            with session.call_value("num_native_qlinalg_det", normalized) as determinant:
                assert int(determinant) == 4 * wide - 2
                print(f"Exact determinant: {determinant.text()}")

    with rssn.Value.tensor_cells(session, [Fraction(1, 3)], ()) as scalar:
        assert scalar.json()["tensor"]["shape"] == []
        assert scalar[0].text() == "1/3"
    with rssn.Value.tensor_cells(session, [], (2, 0, 3), exemplar=0) as empty:
        assert empty.json()["tensor"]["shape"] == [2, 0, 3]
        assert empty.json()["tensor"]["elem"] == "int"

    for values, shape in [([[1], [2, 3]], (2,)), ([1], (-1,)), ([1], (2, 2))]:
        try:
            rssn.Value.tensor_cells(session, values, shape)
        except ValueError:
            pass
        else:
            raise AssertionError("malformed tensor input was accepted")

    with rssn.Session() as foreign, rssn.Value.of(foreign, wide) as cell:
        try:
            rssn.Value.tensor_cells(session, [], (0, 1), exemplar=cell)
        except ValueError as error:
            assert "another session" in str(error)
        else:
            raise AssertionError("a foreign tensor cell was accepted")
