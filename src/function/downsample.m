function f = downsample(f_big, os)
  % f_big  big image
  % os     oversampling factor

  f = zeros(size(f_big)/os);
  os2 = ceil(os/2);
  
  switch (ndims(f_big))
    case 2
      for i1 = 1:2*os2-1
        for i2 = 1:2*os2-1
          f(2:end,2:end) = f(2:end,2:end) + ...
            f_big(os2+i1:os:end-os+os2, os2+i2:os:end-os+os2)/os^2;
        end
      end
      
    case 3
      os2 = ceil(os/2);
      for i1 = 1:2*os2-1
        for i2 = 1:2*os2-1
          for i3 = 1:2*os2-1
            f(2:end,2:end,2:end) = f(2:end,2:end,2:end) + ...
              f_big(os2+i1:os:end-os+os2, os2+i2:os:end-os+os2, os2+i3:os:end-os+os2)/os^3;
          end
        end
      end
      % for i1 = 0:os-1
      %   for i2 = 0:os-1
      %     for i3 = 0:os-1
      %       f_smooth(2:end,2:end,2:end) = f_smooth(2:end,2:end,2:end) + ...
      %         f_big(os+i1:os:end-1, os+i2:os:end-1, os+i3:os:end-1)/os^3;
      %     end
      %   end
      % end
    otherwise
      error('f must be two- or three-dimensional')
  end
end
