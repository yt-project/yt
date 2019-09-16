from yt.testing import fake_random_ds
import yt.utilities.answer_testing.framwork as fw
from yt.utilities.answer_testing import utils


# Answer file
answer_file = 'connected_sets.yaml'


class TesetConnectedSets(fw.AnswerTest):
    def test_connected_sets(self):
        ds = fake_random_ds(16, nprocs=8, particles=16**3)
        data_source = ds.disk([0.5, 0.5, 0.5], [0., 0., 1.], (8, 'kpc'), (1, 'kpc'))
        field = ("gas", "density")
        min_val, max_val = data_source[field].min()/2, data_source[field].max()/2
        data_source.extract_connected_sets(field, 5, min_val, max_val,
                                           log_space=True, cumulative=True)
        hd = OrderedDict()
        ecs_hd = utils.generate_hash(
            self.extract_connected_sets_test(ds, data_source, field, 5, min_val, max_val)
        )
        hd['extract_connected_sets'] = ecs_hd
        hd = {'test_connected_sets' : hd}
        utils.handle_hashes(self.save_dir, answer_file, hd, self.answer_store)
