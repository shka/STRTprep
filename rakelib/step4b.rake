##
## Step 4b - Count in the regions
##

def step4b_sources(path)
  return ['out/byTFE/regions.bed.gz',
          path.sub(/^tmp\/byTFE/, 'tmp').sub(/\.step4b$/, '.step2e')]
end

rule /^tmp\/byTFE\/[^\/]+\.step4b$/ => [->(path){ step4b_sources(path) }] do |t|
  step3b_job(t)
end

#

rule '.step4b_cnt' => '.step4b' do |t|
  step3b_cnt_job(t)
end

#

task :clean_step4b do
  LIBIDS.each do |libid|
    sh "rm -rf tmp/byTFE/#{libid}.*.step4b"
    sh "rm -rf tmp/byTFE/#{libid}.*.step4b_cnt"
  end
end
