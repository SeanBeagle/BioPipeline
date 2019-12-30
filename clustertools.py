import os
import subprocess


class Module:
    """Environment Modules
    Load environmental modules in linux shell.
    """
    CMD = "/usr/local/Modules/default/bin/modulecmd"

    @classmethod
    def load(cls, *modules):
        """Loads environmental modules."""
        for module in modules:
            try:
                exec(subprocess.Popen([cls.CMD, 'python', 'load', module], stdout=subprocess.PIPE).communicate()[0])
            except Exception as e:
                print(e)


class LSF:
    """Load Sharing Facility
    Submit distributed jobs to scheduler.
    """
    @staticmethod
    def bsub(command, shell=False):
        """Submit job on HPC."""
        if isinstance(command, str):
            command = command.split()
        if isinstance(command, list):
            if shell:
                command = f"bsub '{' '.join(command)}'"
                process = subprocess.Popen(command, shell=True)
            else:
                command = ['bsub'] + command
                process = subprocess.run(command, capture_output=True)
        else:
            return
        # id = str(process.stdout).split("<")[1].split(">")[0].strip()
        # subprocess.run(f"bjobs job_ID {id}")
        # return {'pid': id, 'command': command}

    @staticmethod
    def run(command, shell=False):
        if isinstance(command, str):
            command = command.split()
        if isinstance(command, list):
            if shell:
                command = ' '.join(command)
                return subprocess.Popen(command, shell=True)
            else:
                command = command
                return subprocess.run(command, capture_output=True)

