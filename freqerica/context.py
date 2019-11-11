import subprocess

cmd = 'pip show qulacs | grep Version:'
result = (subprocess.Popen(cmd, stdout=subprocess.PIPE,
                           shell=True).communicate()[0]).decode('utf-8')

# result == "Version: X.X.X.X\n"
version = result[9:-1]

rotation_factor = 2.

if version=='0.0.4':
    rotation_factor = 1.

print('Using qulacs', version, '   rotation_factor =', rotation_factor)
