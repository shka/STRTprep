require 'yaml'

PROCS = `gnproc`.to_i

tmp = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'conf.yaml')
CONF = tmp['LIBS']

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

# targets = Array.new
# LIBIDS.each do |libid|
#   WELLS.each do |wellid|
#     targets.push("tmp/cg/#{libid}.#{wellid}.step3c")
#   end
# end

task :v2 => ['out/cg/nreads.txt.gz']

task :default => :v2
