from unittest import TestCase
from unittest.mock import Mock

from make_prg.update.update_shared_data import SingletonUpdateSharedData


class SingletonUpdateSharedDataTest(TestCase):
    def test___SingletonUpdateSharedData_is_created_and_init_only_once(self):
        denovo_variants_db = Mock()
        aligner = Mock()

        # creates and init the singleton
        update_shared_data = SingletonUpdateSharedData(denovo_variants_db, aligner)

        # checks if object is initialised correctly
        self.assertEqual(update_shared_data.data.denovo_variants_db, denovo_variants_db)
        self.assertEqual(update_shared_data.data.aligner, aligner)

        # gets the singleton again
        update_shared_data = SingletonUpdateSharedData()

        # checks that object is not re-initialised, and still contains the same contents
        self.assertEqual(update_shared_data.data.denovo_variants_db, denovo_variants_db)
        self.assertEqual(update_shared_data.data.aligner, aligner)

        # gets the singleton again
        update_shared_data = SingletonUpdateSharedData()

        # checks that object is not re-initialised, and still contains the same contents
        self.assertEqual(update_shared_data.data.denovo_variants_db, denovo_variants_db)
        self.assertEqual(update_shared_data.data.aligner, aligner)
