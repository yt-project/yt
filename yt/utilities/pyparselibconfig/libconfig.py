from __future__ import print_function

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

class libconfig(dict):
    def __init__(self, config=None):
        if config is not None:
            self.read(config)

    def read(self, config):
        if not hasattr(config, "read"):
            cfile = open(config, 'r')
        else:
            cfile = config

        # Strip out spaces and blanks
        lines = [line.strip() for line in cfile if len(line) > 0]

        # Strip out comments
        lines = [line for line in lines if not (line.startswith('#') or
                                                line.startswith('/'))]

        # Concatenate
        oneline = ''
        for line in lines:
            oneline += line

        statements = oneline.split(';')

        self.parse_statements(self, statements)

    def parse_statements(self, this_dict, statements):
        while len(statements) > 0:
            statement = statements.pop(0)
            if len(statement) == 0:
                continue
            if ':' in statement:
                # DICTIONARY
                new_level_lines = []
                k, v = statement.split(':', 1)
                k = k.strip(' :')
                v = v.strip(' {')
                level = 1 + v.count(':')
                this_dict[k] = {}
                new_level_lines.append(v)

                while(level > 0 and len(statements) > 0):
                    nextline = statements.pop(0).strip()
                    level += nextline.count('{')
                    if nextline == '}' and level == 1:
                        level = 0
                        break
                    new_level_lines.append(nextline)
                    level -= nextline.count('}')
                self.parse_statements(this_dict[k], new_level_lines)
            else:
                k, v = statement.split('=')
                k = k.strip()
                v = v.strip()
                this_dict[k] = self.correct_type(v)

    def correct_type(self, v):
        if v == "true":
            v = "True"
        elif v == "false":
            v = "False"
        # ...are we really evaling this?  We should work around that somehow.
        return eval(v)

    def write(self, filename):
        f = open(filename, 'w')

        self.write_dict(f, self, 0)

    def write_dict(self, f, this_dict, level):
        tab = ' '*4

        dict_dict = {}
        for k, v in this_dict.items():
            if type(v) == dict:
                dict_dict[k] = v
            else:
                if type(v) == str:
                    v = '"%s"' % v
                f.writelines(tab*level + '%s = %s;\n' % (k, v))

        for k, v in dict_dict.items():
            f.writelines('\n')
            f.writelines(tab*level + '%s :\n' % k)
            f.writelines(tab*level+'{\n')
            self.write_dict(f, v, level+1)
            f.writelines(tab*level+'};\n')

if __name__ == '__main__':
    cfg = libconfig()
    cfg.read('test_config.cfg')
    print(cfg)
