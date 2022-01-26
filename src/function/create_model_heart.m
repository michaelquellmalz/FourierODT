function [f cp] = create_model_heart(x,y,a,c0,fscale)
% a=8;
% [x,y]=meshgrid(linspace(-10,10,200));
  f = zeros(size(x));
  f(x.^2/a^2 + y.^2/a^2 <= 1) = 0.25*fscale; % 4*
  
  % Heart
  xx = 4*x/a; yy=4*y/a+1.8;
  f((xx.^2+yy.^2-1).^3 <= xx.^2.*yy.^3) = 0.50*fscale; % 4*
  
  % Ellipse
  xx = 5*(x)/a; yy=5*(y)/a-1.8;
  f(xx.^2 + yy.^2 <= 1) = 0.50;%*fscale;
  
  % convert to wave speed: 
  freq=1;
  omega = 2*pi*freq ; 
  cp = omega ./ sqrt((f+omega^2/c0^2)) ; clear omega
%   fprintf('********************** at frequency %f Hz for c0=%f\n', freq,c0);
%   fprintf('********************** min/max for f = [%f %f]\n', min(min(f)),max(max(f)));
%   fprintf('********************** min/max for c = [%f %f]\n', min(min(cp)),max(max(cp)));
%   figure; 
%   subplot(1,2,1); imagesc(real(f))  ; axis image; colormap jet(256) ; colorbar ; set(gca,'ydir','reverse') ; hold on ; title('f');
%   subplot(1,2,2); imagesc(real(cp)) ; axis image; colormap jet(256) ; colorbar ; set(gca,'ydir','reverse') ; hold on ; title('cp');
  
end

% surf(x,y,f,'EdgeColor','none')
% colorbar
% axis equal
% view(2)
