##
## Step 3c - Count in the regions
##

def step3c_sources(path)
  return ['tmp/cg/regions_forQC.bed.gz',
          path.sub(/^tmp\/cg/, 'tmp').sub(/\.step3c$/, '.step2e')]
end

def step3c_job(t)
  step3b_job(t)
end

rule /^tmp\/cg\/[^\/]+\.step3c$/ => [->(path){ step3c_sources(path) }] do |t|
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
    sh "rm -rf tmp/cg/#{libid}.*.step3c"
    sh "rm -rf tmp/cg/#{libid}.*.step3c_cnt"
  end
end
