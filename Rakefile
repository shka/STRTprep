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

task :qc => qc_targets + ['out/byGene/samples.xls']

##

classes = ['global']
begin
  infp = open('src/samples.csv', 'rt')
  colnames = infp.gets.rstrip.split(',')
  infp.close
  (colnames.select { |colname| /^CLASS\.\d+$/.match(colname) }).each do |cls|
    classes.push(/^CLASS\.(\d+)$/.match(cls).to_a[1])
  end
end

plugin_byGene_targets = Array.new
plugin_byTFE_targets = Array.new

Dir.glob('plugins/*') do |script|
  if /[\#~]$/ !~ script
    plugin = script.pathmap('%n')
    classes.each do |cls|
      target = "out/byGene/plugin_#{plugin}_#{cls}.timestamp"
      plugin_byGene_targets.push(target)
      file target => ['out/byGene/diffexp.csv', 'out/byGene/samples.csv'] do |t|
        sh "#{script} byGene #{cls} #{t.sources.join(' ')} conf.yaml > #{t.name.pathmap("%X")}.log 2>&1"
        sh "touch #{t.name}"
      end
      target = "out/byTFE/plugin_#{plugin}_#{cls}.timestamp"
      plugin_byTFE_targets.push(target)
      file target => ['out/byTFE/diffexp.csv', 'out/byGene/samples.csv'] do |t|
        sh "#{script} byTFE #{cls} #{t.sources.join(' ')} conf.yaml > #{t.name.pathmap("%X")}.log 2>&1"
        sh "touch #{t.name}"
      end
    end
  end
end

task :gene => qc_targets +
              plugin_byGene_targets +
              ['out/byGene/samples.xls',
               'out/byGene/diffexp.xls',
               'out/web/regions_byGene.bed.gz']

##

task :default => qc_targets +
                 plugin_byGene_targets +
                 plugin_byTFE_targets +
                 ['out/byGene/samples.xls',
                  'out/byGene/diffexp.xls',
                  'out/byTFE/diffexp.xls',
                  :web]
