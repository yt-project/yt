import numpy as np

from yt.funcs import mylog

from .parallel_analysis_interface import (
    ResultsStorage,
    _get_comm,
    communication_system,
    parallel_capable,
)

messages = dict(
    task=dict(msg="next"),
    result=dict(msg="result"),
    task_req=dict(msg="task_req"),
    end=dict(msg="no_more_tasks"),
)


class TaskQueueNonRoot:
    def __init__(self, tasks, comm, subcomm):
        self.tasks = tasks
        self.results = {}
        self.comm = comm
        self.subcomm = subcomm

    def send_result(self, result):
        new_msg = messages["result"].copy()
        new_msg["value"] = result
        if self.subcomm.rank == 0:
            self.comm.comm.send(new_msg, dest=0, tag=1)
        self.subcomm.barrier()

    def __next__(self):
        msg = messages["task_req"].copy()
        if self.subcomm.rank == 0:
            self.comm.comm.send(msg, dest=0, tag=1)
            msg = self.comm.comm.recv(source=0, tag=2)
        msg = self.subcomm.bcast(msg, root=0)
        if msg["msg"] == messages["end"]["msg"]:
            mylog.debug("Notified to end")
            raise StopIteration
        return msg["value"]

    # For Python 2 compatibility
    next = __next__

    def __iter__(self):
        return self

    def run(self, callable):
        for task in self:
            self.send_result(callable(task))
        return self.finalize()

    def finalize(self, vals=None):
        return self.comm.comm.bcast(vals, root=0)


class TaskQueueRoot(TaskQueueNonRoot):
    def __init__(self, tasks, comm, njobs):
        self.njobs = njobs
        self.tasks = tasks
        self.results = {}
        self.assignments = {}
        self._notified = 0
        self._current = 0
        self._remaining = len(self.tasks)
        self.comm = comm
        # Set up threading here
        # self.dist = threading.Thread(target=self.handle_assignments)
        # self.dist.daemon = True
        # self.dist.start()

    def run(self, func=None):
        self.comm.probe_loop(1, self.handle_assignment)
        return self.finalize(self.results)

    def insert_result(self, source_id, rstore):
        task_id = rstore.result_id
        self.results[task_id] = rstore.result

    def assign_task(self, source_id):
        if self._remaining == 0:
            mylog.debug("Notifying %s to end", source_id)
            msg = messages["end"].copy()
            self._notified += 1
        else:
            msg = messages["task"].copy()
            task_id = self._current
            task = self.tasks[task_id]
            self.assignments[source_id] = task_id
            self._current += 1
            self._remaining -= 1
            msg["value"] = task
        self.comm.comm.send(msg, dest=source_id, tag=2)

    def handle_assignment(self, status):
        msg = self.comm.comm.recv(source=status.source, tag=1)
        if msg["msg"] == messages["result"]["msg"]:
            self.insert_result(status.source, msg["value"])
        elif msg["msg"] == messages["task_req"]["msg"]:
            self.assign_task(status.source)
        else:
            mylog.error("GOT AN UNKNOWN MESSAGE: %s", msg)
            raise RuntimeError
        if self._notified >= self.njobs:
            raise StopIteration


def task_queue(func, tasks, njobs=0):
    comm = _get_comm(())
    if not parallel_capable:
        mylog.error("Cannot create task queue for serial process.")
        raise RuntimeError
    my_size = comm.comm.size
    if njobs <= 0:
        njobs = my_size - 1
    if njobs >= my_size:
        mylog.error(
            "You have asked for %s jobs, but only %s processors are available.",
            njobs,
            (my_size - 1),
        )
        raise RuntimeError
    my_rank = comm.rank
    all_new_comms = np.array_split(np.arange(1, my_size), njobs)
    all_new_comms.insert(0, np.array([0]))
    for i, comm_set in enumerate(all_new_comms):
        if my_rank in comm_set:
            my_new_id = i
            break
    subcomm = communication_system.push_with_ids(all_new_comms[my_new_id].tolist())

    if comm.comm.rank == 0:
        my_q = TaskQueueRoot(tasks, comm, njobs)
    else:
        my_q = TaskQueueNonRoot(None, comm, subcomm)
    communication_system.pop()
    return my_q.run(func)


def dynamic_parallel_objects(tasks, njobs=0, storage=None, broadcast=True):
    comm = _get_comm(())
    if not parallel_capable:
        mylog.error("Cannot create task queue for serial process.")
        raise RuntimeError
    my_size = comm.comm.size
    if njobs <= 0:
        njobs = my_size - 1
    if njobs >= my_size:
        mylog.error(
            "You have asked for %s jobs, but only %s processors are available.",
            njobs,
            (my_size - 1),
        )
        raise RuntimeError
    my_rank = comm.rank
    all_new_comms = np.array_split(np.arange(1, my_size), njobs)
    all_new_comms.insert(0, np.array([0]))
    for i, comm_set in enumerate(all_new_comms):
        if my_rank in comm_set:
            my_new_id = i
            break
    subcomm = communication_system.push_with_ids(all_new_comms[my_new_id].tolist())

    if comm.comm.rank == 0:
        my_q = TaskQueueRoot(tasks, comm, njobs)
        my_q.comm.probe_loop(1, my_q.handle_assignment)
    else:
        my_q = TaskQueueNonRoot(None, comm, subcomm)
        if storage is None:
            for task in my_q:
                yield task
        else:
            for task in my_q:
                rstore = ResultsStorage()
                yield rstore, task
                my_q.send_result(rstore)

    if storage is not None:
        if broadcast:
            my_results = my_q.comm.comm.bcast(my_q.results, root=0)
        else:
            my_results = my_q.results
        storage.update(my_results)

    communication_system.pop()
