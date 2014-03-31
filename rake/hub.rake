require 'parallel'

####
#
# dummy paths
#
tmp = Array.new
LIBIDS.each { |libid|
  tmp2 = Array.new
  open(File.expand_path(CONF[libid]['LAYOUT'])).each { |line|
    cols = line.rstrip.split
    wellid = cols[0]
    ['fwd', 'rev'].each { |str|
      tmp2.push("out/hub/#{libid}.#{wellid}.#{str}.uniq.5p.bw")
    }
  }
  #
  timestamp = "tmp/#{libid}.hub.timestamp"
  file timestamp => tmp2 do |t|
    sh "touch #{t.name}"
    begin
      Parallel.each(t.prerequisites, :in_threads => PROCS) { |target|
        Rake::Task[target].invoke
      }
    rescue
      sh "rm -rf #{timestamp}"
    end
  end
  tmp.push(timestamp)
}
task :hub => tmp

####
#
# file out/hub/#{libid}.#{wellid}.(fwd|rev).bam
#      => out/ali/#{libid}.#{wellid}/accepted_hits.bam
#
rule /^out\/hub\/[^.]+\.[^.]+\.(fwd|rev)\.bam/ => proc { |target|
  target.sub('hub', 'ali').sub(/.(fwd|rev)/, '/accepted_hits')
} do |t|
  libid, wellid, str = /\/([^.\/]+)\.([^.]+)\.(fwd|rev)\.bam/.match(t.name).to_a[1..3]
  #
  orgbam = t.prerequisites[0]
  strbam = t.name
  struniq5pbed = strbam.sub('out/', 'tmp/').sub('.bam', '.uniq.5p.bed.gz')
  struniq5pwig = struniq5pbed.sub('.bed', '.wig')
  struniq5pbigwig = struniq5pwig.sub('.wig.gz', '.bw').sub('tmp/', 'out/')
  fai = File.expand_path(CONF[libid]['GENOMESPIKERIBO']+'.fa.fai')
  #
  sh "mkdir -p out/hub tmp/hub"
  begin
    sh "samtools view -h #{orgbam} | grep -P '(^@|XS:A:\\#{str == 'fwd' ? '+' : '-'})' | samtools view -S -b - > #{strbam}"
  rescue
    puts "... no alignments for bam on #{libid}.#{wellid}.#{str}; skipped"
    sh "rm -rf #{struniq5pbed}; touch #{struniq5pbed}"
  else
    sh "bin/_process_bam2bed_uniq_5p.rb #{str} #{orgbam} #{strbam} | gzip -c > #{struniq5pbed}"
    sh "bin/_process_bed2wig_5p.rb #{libid}.#{wellid} #{str} #{struniq5pbed} | gzip -c > #{struniq5pwig}"
    begin
      sh "wigToBigWig #{struniq5pwig} #{fai} #{struniq5pbigwig}"
    rescue
      puts "... no alignments for bigWig on #{libid}.#{wellid}.#{str}; skipped"
    end
  end
end

####
#
# file tmp/hub/#{libid}.#{wellid}.(fwd|rev).uniq.bed.gz
#      => [ out/ali/#{libid}.#{wellid}/accepted_hits.bam,
#           out/hub/#{libid}.#{wellid}.(fwd|rev).bam ]
#
rule /^tmp\/hub\/[^.]+\.[^.]+\.(fwd|rev)\.uniq\.5p\.bed\.gz/ => proc { |target|
  [ target.sub('tmp/hub', 'out/ali').sub(/\.(fwd|rev)\.uniq\.5p\.bed\.gz/, '/accepted_hits.bam'),
    target.sub('tmp/', 'out/').sub('.uniq.5p.bed.gz', '.bam') ]
}

####
#
# file tmp/hub/#{libid}.#{wellid}.(fwd|rev).uniq.5p.wig.gz
#      => tmp/hub/#{libid}.#{wellid}.(fwd|rev).uniq.5p.bed.gz
#
rule /^tmp\/hub\/[^.]+\.[^.]+\.(fwd|rev)\.uniq\.5p.wig\.gz/ => proc { |target|
  target.sub('.wig', '.bed')
}

####
#
# file out/hub/#{libid}.#{wellid}.(fwd|rev).uniq.5p.bw
#      => tmp/hub/#{libid}.#{wellid}.(fwd|rev).uniq.5p.wig.gz
#
rule /^out\/hub\/[^.]+\.[^.]+\.(fwd|rev)\.uniq\.5p\.bw/ => proc { |target|
  [ target.sub('out', 'tmp').sub('.bw', '.wig.gz'),
    File.expand_path(CONF[/hub\/([^.]+)\./.match(target).to_a[1]]['GENOMESPIKERIBO']+'.fa.fai') ]
}

####
#
# cleaning
#
task 'clean_hub' do
  LIBIDS.each { |libid| sh "rm -rf tmp/hub/#{libid}.* out/hub/#{libid}.* tmp/#{libid}.hub.timestamp" }
end
