## Simple routines to write parameter files


# ----------------------------------------------------------------------
# functions to write RASCAS parameter files from dictionaries.
# ----------------------------------------------------------------------

class param_section():
    def __init__(self,name,params):
        self.name = name
        self.params = params


def write_param_section(f,s):
    f.write('[%s] \n'%s.name)
    for k in s.params.keys():
        f.write('  %s = %s \n'%(k,s.params[k]))
    f.write('\n') # skip a line between sections
        
def write_parameter_file(filename,params):
    f = open(filename,'w')
    for s in params:
        write_param_section(f,s)
    f.close()
    
def print_param_section(ps):
    print('[%s] \n'%ps.name)
    for k in ps.params.keys():
        print('  %s = %s \n'%(k,ps.params[k]))
        
        
