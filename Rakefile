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

quick_targets = Array.new
LIBIDS.each do |libid|
  CONF[libid]['FASTQS'].each_index do |runid|
    quick_targets.push("tmp/#{libid}.#{runid}.step1b")
  end
end
LIBWELLIDS.each do |libwellid|
  quick_targets.push("tmp/#{libwellid}.step2g")
end
LIBWELLIDS.each do |libwellid|
  quick_targets.push("tmp/byGene/#{libwellid}.step3b",
                     "tmp/byGene/#{libwellid}.step3c")
end

task :quick => quick_targets + ['out/byGene/fluctuation.txt.gz', 'out/byGene/samples.csv']

##

full_targets = Array.new
LIBWELLIDS.each do |libwellid|
  quick_targets.push("tmp/byGene/#{libwellid}.step4b")
end

task :default => quick_targets + full_targets + ['out/byGene/diffexp.csv', 'out/byTFE/diffexp.csv', :web]
