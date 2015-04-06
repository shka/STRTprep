##
## Step 4a - transcript assembl
##

require 'csv'

def step4a_sources(t)
  cls = t.pathmap("%d").pathmap("%f").pathmap("%X")
  src = ['src/samples.csv']
  samples = CSV.table(src[0])
  samples.each do |row|
    if row[:classtfe] == cls
      src.push("out/bam/#{row[:library]}.#{row[:well]}.bam")
    end
  end
  return src
end

rule /\.step4a\/transcripts\.gtf$/ => [->(path){ step4a_sources(path) }] do |t|
  dir = t.name.pathmap("%d")
  mkdir_p dir
  cls = t.name.pathmap("%d").pathmap("%f").pathmap("%X")
  sh "samtools merge -f #{dir}/merged.bam #{t.sources[1..-1].join(' ')} > #{t.name}.log 2>&1"
  sh "(cufflinks -o #{dir} -p #{PROCS} --library-type fr-secondstrand -L #{cls} #{dir}/merged.bam) >> #{t.name}.log 2>&1"
end

##

step4a_targets = Array.new
begin 
  samples = CSV.table('src/samples.csv')
  classes = Array.new
  samples.each { |row| classes.push(row[:classtfe]) if row[:classtfe] != 'NA' }
  classes.uniq.each do |cls|
    step4a_targets.push("tmp/#{cls}.step4a/transcripts.gtf")
  end
end

task :step4a => step4a_targets
