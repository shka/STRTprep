require 'csv'
require 'pp'
require 'yaml'

PROCS = ENV.key?('PROCS') ? ENV['PROCS'].to_i : `gnproc`.to_i
THREADS = Rake.application.thread_pool.statistics[:max_active_threads]+1

conf = YAML.load_file(ENV.key?('CONF') ? ENV['CONF'] : 'src/conf.yaml')
LIBRARIES = conf['LIBRARIES']
PLUGINS = conf['PLUGINS']
PREPROCESS = conf['PREPROCESS']

LIBIDS = ENV.key?('LIBS') ? ENV['LIBS'].split(',') : LIBRARIES.keys

def getLayout(libid)
  return LIBRARIES[libid].key?('LAYOUT') ? LIBRARIES[libid]['LAYOUT'] : PREPROCESS['LAYOUT']
end

LIBWELLIDS = Array.new
LIBIDS.each do |libid|
  layoutFile = getLayout(libid)
  open(layoutFile).each do |line|
    wellid, barcode = line.rstrip.split(/\t/)
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

## QC tasks

### QC1: parallel by lanes

qcTargets1 = Array.new

LIBIDS.each do |libid|
  LIBRARIES[libid]['FASTQS'].each_index do |runid|
    qcTargets1.push("tmp/#{libid}.#{runid}.step1b")
  end
end

task :qc1 => qcTargets1

#### QC2: parallel by libraries

qcTargets2 = Array.new

LIBIDS.each do |libid|
  qcTargets2.push("tmp/#{libid}.step2a", "out/fig.#{libid}.qualifiedReads.pdf")
end

task :qc2 => qcTargets2

#### QC3: parallel by project

qcTargets3 = ['tmp/step2c/accepted_hits.samUniqSortedByAcc']

task :qc3 => qcTargets3

#### QC4: parallel by wells

qcTargets4 = Array.new

LIBWELLIDS.each do |libwellid|
  qcTargets4.push("tmp/#{libwellid}.step2d_cnt",
                  "tmp/#{libwellid}.step2f_cnt",
                  "tmp/#{libwellid}.step2g",
                  "tmp/byGene/#{libwellid}.step3b_cnt",
                  "tmp/byGene/#{libwellid}.step3c_cnt",
                  "out/seq/#{libwellid}.fq.gz")
end

multitask :qc4 => qcTargets4

####

task :qc => [:qc1, :qc2, :qc3, :qc4, 'out/byGene/samples.xls']

## TFE intermediate tasks

#### TFE1: parallel by project

tfeTargets1 = Array.new
begin
  samples = CSV.table('src/samples.csv')
  classes = Array.new
  samples.each { |row| classes.push(row[:classtfe]) if row[:classtfe] != 'NA' }
  classes.uniq.each do |cls|
    tfeTargets1.push("tmp/byTFE/#{cls}.step4a/transcripts.gtf")
  end
end
tfeTargets1 << 'out/byTFE/regions.bed.gz'

task :tfe1 => tfeTargets1

#### TFE2: parallel by wells

tfeTargets2 = Array.new
LIBWELLIDS.each do |libwellid|
  tfeTargets2.push("tmp/byTFE/#{libwellid}.step4b_cnt")
end

multitask :tfe2 => tfeTargets2

##

classes = ['global']
begin
  infp = open('src/samples.csv', 'rt')
  colnames = infp.gets.rstrip.split(',')
  infp.close
  (colnames.select { |colname| /^CLASS\.\d+$/.match(colname) }).each do |cls|
    classes.push(/^CLASS\.(\d+)$/.match(cls).to_a[1].to_i)
  end
end

plugin_byGene_targets = Array.new
plugin_byTFE_targets = Array.new

if !PLUGINS.nil?
  PLUGINS.each_key do |plugin|
    script = "plugins/#{plugin}.R"
    if File.exist?(script)
      classes.each do |cls|
        if PLUGINS[plugin].key?(cls)
          target = "out/byGene/plugin_#{plugin}_#{cls}.timestamp"
          plugin_byGene_targets.push(target)
          file target => ['out/byGene/diffexp.csv',
                          'out/byGene/samples.csv',
                          'src/conf.yaml'] do |t|
            sh "#{script} byGene #{cls} #{t.sources.join(' ')} > #{t.name.pathmap("%X")}.log 2>&1"
            sh "touch #{t.name}"
          end
          target = "out/byTFE/plugin_#{plugin}_#{cls}.timestamp"
          plugin_byTFE_targets.push(target)
          file target => ['out/byTFE/diffexp.csv',
                          'out/byGene/samples.csv',
                          'src/conf.yaml'] do |t|
            sh "#{script} byTFE #{cls} #{t.sources.join(' ')} > #{t.name.pathmap("%X")}.log 2>&1"
            sh "touch #{t.name}"
          end
        end
      end
    end
  end
end

##

task :gene => [:qc,
               'out/byGene/samples.xls',
               'out/byGene/diffexp.xls',
               'out/web/regions_byGene.bed.gz']

task :plugins_gene => [:gene] + plugin_byGene_targets

task :tfe => [:qc,
              'out/byGene/samples.xls',
              'out/byTFE/diffexp.xls']

task :plugins_tfe => [:tfe] + plugin_byTFE_targets

task :default => [:gene, :plugins_gene, :tfe, :plugins_tfe, :web]

task :check do
  pp conf
end
