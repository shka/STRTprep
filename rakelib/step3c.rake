##
## Step 3c - Count in the regions for QC
##

def step3c_sources(path)
  return ['tmp/byGene/regions_forQC.bed.gz',
          path.sub(/^tmp\/byGene/, 'tmp').sub(/\.step3c$/, '.step2e')]
end

def step3c_job(t)
  step3b_job(t)
end

rule /^tmp\/byGene\/[^\/]+\.step3c$/ => [->(path){ step3c_sources(path) }] do |t|
  step3c_job(t)
end

#

def step3c_cnt_job(t)
  step3b_cnt_job(t)
end

rule '.step3c_cnt' => '.step3c' do |t|
  step3c_cnt_job(t)
end

#

task :clean_step3c do
  LIBIDS.each do |libid|
    sh "rm -rf tmp/byGene/#{libid}.*.step3c"
    sh "rm -rf tmp/byGene/#{libid}.*.step3c_cnt"
  end
end
