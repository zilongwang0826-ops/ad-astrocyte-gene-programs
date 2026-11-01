"""Tiny fixtures exercise raw storage, donor metadata and count conservation."""
from contextlib import contextmanager
from pathlib import Path
import sys
import tempfile
import unittest

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from data_utils import align_barcodes, aggregate_block, blocks, load_metadata


@contextmanager
def stored_csc(data, indices, indptr, shape):
    """Use an actual HDF5 group so slicing and on-disk dtypes are exercised."""
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'counts.h5'
        with h5py.File(path, 'w') as handle:
            group = handle.create_group('matrix')
            for name, value in dict(data=data, indices=indices,
                                    indptr=indptr, shape=shape).items():
                group.create_dataset(name, data=value)
        with h5py.File(path, 'r') as handle:
            yield handle['matrix']


class DataChecks(unittest.TestCase):
    def test_reordered_barcodes(self):
        np.testing.assert_array_equal(
            align_barcodes(['b', 'c', 'a'], pd.DataFrame({'Barcode': ['a', 'b']})),
            [2, 0])

    def test_missing_cell_rejected(self):
        with self.assertRaisesRegex(ValueError, 'absent from expression matrix'):
            align_barcodes(['a'], pd.DataFrame({'Barcode': ['absent']}))

    def test_duplicate_barcode_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicated matrix barcodes'):
            align_barcodes(['a', 'a'], pd.DataFrame({'Barcode': ['a']}))

    def test_donor_sum_conservation(self):
        matrix = sparse.csc_matrix([[1, 2, 100, 3], [2, 0, 100, 5]])
        out = aggregate_block(matrix, ['D2', 'D1', '', 'D2'], ['D1', 'D2'])
        np.testing.assert_array_equal(out, [[2, 4], [0, 7]])
        self.assertEqual(out.sum(), matrix[:, [0, 1, 3]].sum())

    def test_no_selection(self):
        out = aggregate_block(sparse.eye(2, format='csc'), ['', ''], ['D1'])
        np.testing.assert_array_equal(out, [[0], [0]])


class StorageChecks(unittest.TestCase):
    def test_unsorted_storage_and_empty_column_conserve_counts(self):
        expected = np.array([[2, 0, 0, 4], [0, 0, 5, 0], [3, 0, 0, 0]])
        with stored_csc(np.array([3, 2, 5, 4], dtype=np.int32),
                        [2, 0, 1, 0], [0, 2, 2, 3, 4], [3, 4]) as group:
            for size in (1, 2, 3, 8):
                with self.subTest(block_size=size):
                    chunks = list(blocks(group, block_size=size))
                    self.assertEqual(chunks[0][0], 0)
                    self.assertEqual(chunks[-1][1], 4)
                    self.assertTrue(all(m.has_canonical_format for _, _, m in chunks))
                    reconstructed = sparse.hstack([m for _, _, m in chunks]).toarray()
                    np.testing.assert_array_equal(reconstructed, expected)
                    self.assertEqual(int(reconstructed.sum()), 14)
                    donor_counts = sum(
                        (aggregate_block(m, ['D2', 'D1', 'D1', 'D2'][start:stop],
                                         ['D1', 'D2']) for start, stop, m in chunks),
                        np.zeros((3, 2), dtype=np.int64))
                    np.testing.assert_array_equal(donor_counts, [[0, 6], [5, 0], [0, 3]])

    def test_empty_matrix_storage(self):
        with stored_csc(np.array([], dtype=np.int32), np.array([], dtype=np.int64),
                        [0, 0, 0], [3, 2]) as group:
            observed = list(blocks(group, block_size=1))
            self.assertEqual(len(observed), 2)
            self.assertTrue(all(m.nnz == 0 and m.shape == (3, 1) for _, _, m in observed))

    def test_duplicate_coordinates_rejected(self):
        with stored_csc([2, 3], [1, 1], [0, 2], [3, 1]) as group:
            with self.assertRaisesRegex(ValueError, 'Duplicated entries'):
                list(blocks(group))

    def test_malformed_pointers_rejected(self):
        for pointer in ([1, 1, 2], [0, 2, 1], [0, 2], [0, 1, 3]):
            with self.subTest(indptr=pointer):
                with stored_csc([2, 3], [0, 1], pointer, [2, 2]) as group:
                    with self.assertRaisesRegex(ValueError, 'CSC'):
                        list(blocks(group))

    def test_storage_lengths_rejected(self):
        with stored_csc([2, 3], [0], [0, 2], [2, 1]) as group:
            with self.assertRaisesRegex(ValueError, 'storage length mismatch'):
                list(blocks(group))

    def test_invalid_gene_index_rejected(self):
        for index in (-1, 3):
            with self.subTest(index=index):
                with stored_csc([2], [index], [0, 1], [3, 1]) as group:
                    with self.assertRaisesRegex(ValueError, 'Gene index out of bounds'):
                        list(blocks(group))

    def test_negative_and_noninteger_count_storage_rejected(self):
        for count in (np.array([-1], dtype=np.int32), np.array([0.5]), np.array([2.0])):
            with self.subTest(dtype=str(count.dtype), value=count.tolist()):
                with stored_csc(count, [0], [0, 1], [2, 1]) as group:
                    with self.assertRaisesRegex(ValueError, 'nonnegative raw integers'):
                        list(blocks(group))

    def test_noninteger_coordinate_storage_rejected(self):
        cases = [
            ([0.5], [0, 1], [2, 1], 'integer storage'),
            ([0], [0, 0.5, 1], [2, 2], 'CSC pointers'),
            ([0], [0, 1], [2.5, 1], 'integer CSC shape'),
        ]
        for indices, pointers, shape, message in cases:
            with self.subTest(indices=indices, pointers=pointers, shape=shape):
                with stored_csc([3], indices, pointers, shape) as group:
                    with self.assertRaisesRegex(ValueError, message):
                        list(blocks(group))

    def test_unsigned_descending_pointers_rejected(self):
        # np.diff on uint64 wraps a negative difference to a large positive value.
        with stored_csc([2, 3], [0, 1], np.array([0, 2, 1, 2], dtype=np.uint64),
                        [2, 3]) as group:
            with self.assertRaisesRegex(ValueError, 'Invalid CSC pointers'):
                list(blocks(group))


class MetadataChecks(unittest.TestCase):
    @staticmethod
    def rows():
        return pd.DataFrame({
            'Barcode': ['a', 'b', 'c'], 'SampleID': ['D1', 'D1', 'D2'],
            'Diagnosis': ['AD', 'AD', 'Control'], 'Batch': [1, 1, 2],
            'Cell.Type': ['ASC', 'ASC', 'ASC'], 'Age': [90, 90, 80],
            'Sex': ['F', 'F', 'M'], 'PMI': ['NA', '', '2.5'],
            'RIN': [8.1, 8.1, 7.2], 'Plaque.Stage': ['None', 'None', 'Stage C']})

    def read_fixture(self, frame):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'metadata.csv.gz'
            frame.to_csv(path, index=False)
            return load_metadata(path)

    def test_literal_none_preserved_and_missing_values_explicit(self):
        observed = self.read_fixture(self.rows())
        self.assertEqual(observed['Plaque.Stage'].tolist(), ['None', 'None', 'Stage C'])
        self.assertEqual(observed['PMI'].isna().tolist(), [True, True, False])
        self.assertEqual(observed.loc[2, 'PMI'], 2.5)
        self.assertEqual(observed['SampleID'].tolist(), ['D1', 'D1', 'D2'])

    def test_required_field_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Required metadata fields absent'):
            self.read_fixture(self.rows().drop(columns='SampleID'))

    def test_missing_or_duplicated_barcode_rejected(self):
        for barcode in ('', 'a'):
            with self.subTest(barcode=barcode):
                frame = self.rows()
                frame.loc[1, 'Barcode'] = barcode
                with self.assertRaisesRegex(ValueError, 'barcodes missing or duplicated'):
                    self.read_fixture(frame)

    def test_ambiguous_diagnosis_rejected(self):
        frame = self.rows()
        frame.loc[2, 'Diagnosis'] = 'Unknown'
        with self.assertRaisesRegex(ValueError, 'Ambiguous diagnosis'):
            self.read_fixture(frame)

    def test_inconsistent_donor_fields_rejected(self):
        replacements = {'Diagnosis': 'Control', 'Batch': 2, 'Age': 89,
                        'Sex': 'M', 'PMI': '2.0', 'RIN': 7.5}
        for field, value in replacements.items():
            with self.subTest(field=field):
                frame = self.rows()
                frame.loc[1, field] = value
                with self.assertRaisesRegex(ValueError, 'Inconsistent donor metadata: ' + field):
                    self.read_fixture(frame)

    def test_missing_donor_identifier_rejected(self):
        for donor in ('', 'NA', '   '):
            with self.subTest(donor=donor):
                frame = self.rows()
                frame.loc[1, 'SampleID'] = donor
                with self.assertRaisesRegex(ValueError, 'Donor identifiers missing or blank'):
                    self.read_fixture(frame)


if __name__ == '__main__':
    unittest.main()
