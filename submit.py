'''
Parse a json file with job submission information. Submit jobs.
'''

import json, subprocess, os

class Submitter:
    @staticmethod
    def submit_jobs(jobs):
        '''read a json file with job info and submit jobs'''

        # look for a global config file
        global_options = json.load(os.path.join(os.path.dirname(__file__), 'user.cfg'))["options"]

        config = json.load(jobs)
        options = dict(global_options)
        options.update(config["options"])

        if "jobs" not in config:
            raise RuntimeError("no jobs found in jobs file")

        for job in config["jobs"]:
            if "options" in config:
                if job in config["options"]:
                    pass
            

    @classmethod
    def assert_exists(cls, fns):
        '''assert that files with these names exists'''

        # run ls to update existence
        subprocess.check_output(['ls', '-lah'])

        # check if the files exist
        missing_files = " ".join([fn for fn in fns if not os.path.isfile(fn)])
        if missing_files != "":
            raise RuntimeError("file(s) missing: %s" % missing_files)

    @classmethod
    def assert_doesnt_exist(cls, fns):
        '''assert the files with these names do not exist'''

        # run ls to update existence
        subprocess.check_output(['ls', '-lah'])

        # check if the files exist
        colliding_files = " ".join([fn for fn in fns if os.path.isfile(fn)])
        if colliding_files != "":
            raise RuntimeError("file collision(s): %s" % colliding_files)        

    @classmethod
    def assert_nonempty(cls, fns):
        '''assert that files with these names are nonempty'''

        # check if these files exist first
        cls.assert_exists(fns)

        # check if files are nonempty
        empty_files = " ".join([fn for fn in fns if os.stat(fn).st_size == 0])
        if empty_files != "":
            raise RuntimeError("file(s) empty: %s" % empty_files)
