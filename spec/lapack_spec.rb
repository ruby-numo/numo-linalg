require 'numo/narray'
require File.join(File.dirname(__FILE__), "../ext/numo/lapack/lapack")
#Numo::NArray.debug = true

RSpec.configure do |config|
  config.filter_run :focus
  config.run_all_when_everything_filtered = true
end
#context :focus=>true do ... end

describe 'something' do
  it{expect(Numo::LAPACK).to be_a(Module)}
end
