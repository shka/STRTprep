task :clean do 
  LIBIDS.each { |libid|
    sh "rm -rf out/ali/#{libid}.* out/seq/#{libid}.* out/stat/#{libid}.* tmp/ali/#{libid}.* tmp/seq/#{libid}.*"
  }
end
