from subprocess import PIPE, Popen


def run_cmd(cmd, verbose=False, output=False,error=False):
    stream=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = stream.communicate()
    
    stdout=stdout.decode('utf-8')
    stderr=stderr.decode('utf-8')
    
    print(stderr, flush=True)
    
    if verbose:
        print(stdout, flush=True)
        
        
    if output:
        return stdout
    if error:
        return stderr
    