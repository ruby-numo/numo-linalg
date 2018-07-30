# frozen_string_literal: true

require 'spec_helper'

RSpec.describe Numo::Linalg::Autoloader do
  it 'succeses loading backend libraries' do
    expect { described_class.load_library }.to_not raise_error
    expect(described_class.libs).to_not be_nil
  end

  describe 'private method' do
    let(:autoloader) do
      class DummyLoader
        include Numo::Linalg::Autoloader
      end
      DummyLoader.new
    end

    it 'finds BLAS/LAPACK libraries' do
      lapack_libs = autoloader.send(:find_lapack_libs, ['/usr/lib', '/usr/local/lib', '/usr/local/opt/lapack/lib'])
      expect(lapack_libs[:blas]).to_not be_nil
      expect(lapack_libs[:cblas]).to_not be_nil
      expect(lapack_libs[:lapack]).to_not be_nil
      expect(lapack_libs[:lapacke]).to_not be_nil
    end
  end
end
