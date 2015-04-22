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

require 'csv'
require 'spreadsheet'

rule '.xls' => '.csv' do |t|
  book = Spreadsheet::Workbook.new
  sheet = book.create_worksheet
  
  csv = CSV.read(t.source)
  sheet.row(0).replace csv.first
  
  table = CSV.table(t.source)
  rowidx = 1
  table.each do |row|
    sheet.row(rowidx).replace row.fields
    rowidx += 1
  end
  
  book.write(t.name)
end

##

qc_targets = Array.new
LIBIDS.each do |libid|
  CONF[libid]['FASTQS'].each_index do |runid|
    qc_targets.push("tmp/#{libid}.#{runid}.step1b")
  end
end
LIBIDS.each do |libid|
  qc_targets.push("tmp/#{libid}.step2a")
end
LIBWELLIDS.each do |libwellid|
  qc_targets.push("tmp/#{libwellid}.step2f")
end
LIBWELLIDS.each do |libwellid|
  qc_targets.push("tmp/byGene/#{libwellid}.step3b")
end
LIBWELLIDS.each do |libwellid|
  qc_targets.push("tmp/byGene/#{libwellid}.step3c")
end

task :qc => qc_targets + ['out/byGene/fluctuation.txt.gz',
                          'out/byGene/samples.xls']

##

task :gene => qc_targets + ['out/byGene/diffexp.xls',
                            'out/byGene/samples.xls',
                            'out/web/regions_byGene.bed.gz']

##

full_targets = Array.new
LIBWELLIDS.each do |libwellid|
  full_targets.push("tmp/byTFE/#{libwellid}.step4b")
end

task :default => qc_targets + full_targets + ['out/byGene/diffexp.xls',
                                              'out/byGene/samples.xls',
                                              'out/byTFE/diffexp.xls', :web]
