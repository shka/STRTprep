require 'yaml'

PROCS = `gnproc`.to_i

tmp = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'conf.yaml')
CONF = tmp['LIBS']
DEFAULTS = tmp['DEFAULTS']

LIBIDS = ENV.key?('LIBS') ? ENV['LIBS'].split(',') : CONF.keys

WELLIDS = Array.new
['A', 'B', 'C', 'D', 'E', 'F'].each do |col|
  1.upto(8).each do |row|
    WELLIDS.push("#{col}#{row}")
  end
end

LIBWELLIDS = Array.new
LIBIDS.each do |libid|
  WELLIDS.each do |wellid|
    LIBWELLIDS.push("#{libid}.#{wellid}")
  end
end

##

targets = ['out/cg/diffexp.csv']
LIBWELLIDS.each do |libwellid|
  targets.unshift("tmp/cg/#{libwellid}.step3c_cnt")
end

task :default => targets
