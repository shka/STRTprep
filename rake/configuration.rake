require 'parallel'
require 'yaml'

##

PROCS = Parallel.processor_count

tmp = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'conf.yaml')
CONF = tmp['LIBS']

LIBIDS = ENV.key?('LIBS') ? ENV['LIBS'].split(',') : CONF.keys
