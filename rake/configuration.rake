require 'parallel'
require 'yaml'

##

PROCS = Parallel.processor_count

CONF = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'conf.yaml')

LIBIDS = ENV.key?('LIBS') ? ENV['LIBS'].split(',') : CONF['LIBS'].keys
