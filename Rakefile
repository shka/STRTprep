require 'yaml'

PROCS = `gnproc`.to_i

tmp = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'conf.yaml')
CONF = tmp['LIBS']

LIBIDS = ENV.key?('LIBS') ? ENV['LIBS'].split(',') : CONF.keys

##

targets = Array.new
LIBIDS.each do |libid|
  ['A', 'B', 'C', 'D', 'E', 'F'].each do |col|
    1.upto(8).each do |row|
      targets.push("tmp/cg/#{libid}.#{col}#{row}.step3c")
    end
  end
end

task :v2 => targets

task :default => :v2
