function [f cp] = create_model_shapes(lmax,dx,c0,fscale)
  % clear all;
  % global domain ========================================================
  % dx=0.01; % lmax=5;
  x=[-lmax:dx:lmax];
  nx=length(x); 
  
  f = zeros(nx) ; 
  fmax=0.50*fscale;
  fmid=0.25*fscale;
  % first object is an ellipse.
  ellipse1_x      = 0.0 ; 
  ellipse1_z      = 0.0 ; 
  ellipse1_width  = 4.50; 
  ellipse1_height = 3.25; 
  ellipse1_val    = fmax;
  
  ellipse2in_x     (:) =[ 0.00   1.50   1.50-3.0];  ellipse2out_x     (:) =[ 0.0    1.50     1.50-3.0]; 
  ellipse2in_z     (:) =[-3.00   0.20   0.20+0.2];  ellipse2out_z     (:) =[-4.0  -0.15    -0.15+0.2];
  ellipse2in_width (:) =[ 1.20   0.95   0.95    ];  ellipse2out_width (:) =[ 1.25   0.65     0.65];
  ellipse2in_height(:) =[ 0.90   0.95   0.95    ];  ellipse2out_height(:) =[ 0.75   0.65     0.65];
  ellipse2in_val   (:) =[ fmid   fmid   fmid    ];  
  % two half circles
  halfc_x = [0.00] ; 
  halfc_z = [2.50] ; 
  halfc_r = [1.25] ; 
  halfc_v = [fmid] ; 
  
  
  % ======================================================================
  ixmin=floor((-ellipse1_height+lmax)/dx) + 1 ; 
  ilmax=floor(( ellipse1_height+lmax)/dx) + 1 ; 
  izmin=floor((-ellipse1_width+lmax )/dx) + 1 ; 
  izmax=floor(( ellipse1_width+lmax )/dx) + 1 ; 
  for ix=ixmin:ilmax ; for iz=izmin:izmax
    xx=x(ix) ; zz=x(iz) ; 
  
    cond = (xx-ellipse1_x).^2./ellipse1_height.^2 +  (zz-ellipse1_z).^2 ./ ellipse1_width.^2;
    % allow for some randomness ??
    % tol = -0.20 + (0.20-(-0.20)).*rand(1,1) ; 
    tol = 0.0 ; 
    if(cond <= 1+tol)
      f(ix,iz) = ellipse1_val;  
      
      % Small ellipses inside =============================================
      for k=1:length(ellipse2in_x) 
        ex=ellipse2in_x(k); 
        eh=ellipse2in_height(k);
        ez=ellipse2in_z(k);
        ew=ellipse2in_width(k);
        ev=ellipse2in_val(k) ;
        cond = (xx-ex).^2./eh.^2 +  (zz-ez).^2 ./ ew.^2;
        if(cond <= 1)
          f(ix,iz) = ev;  
        end  
        cond = (xx-ellipse2out_x(k)).^2./ellipse2out_height(k).^2 +  (zz-ellipse2out_z(k)).^2 ./ ellipse2out_width(k).^2;
        if(cond <= 1)
          f(ix,iz) = ellipse1_val; % ellipse2out_val;  
        end  
      end
      % ====================================================================
      
      % half-circles
      for k=1:length(halfc_x)
        cond = sqrt((xx-halfc_x(k))^2 + (zz-halfc_z(k))^2) ;
        if(cond < halfc_r(k)) 
          if(zz-halfc_z(k) >= 0)
            f(ix,iz) = halfc_v(k) ;
          end
        end
      end 
  
    end
  end ; end
  clear xx zz ix iz tol
  % ======================================================================
  
  % convert to wave speed
  % k² = k0² + f 
  %=> (omega/cp) = sqrt(omega² + f)
  %=>  cp = omega ./ sqrt(f+omega²)
  % c0  =1; % 1./(1.+1./3.);
  freq=1;
  omega = 2*pi*freq ; 
  cp = omega ./ sqrt((f+omega^2/c0^2)) ; clear omega
  
%  fprintf('********************** at frequency %f Hz for c0=%f\n', freq,c0);
%   fprintf('********************** min/max for f = [%f %f]\n', min(min(f)),max(max(f)));
%   fprintf('********************** min/max for c = [%f %f]\n', min(min(cp)),max(max(cp)));
       
%   figure; 
%   subplot(1,2,1); imagesc(x,x,real(f)) ; axis image; colormap jet(256) ; colorbar ; set(gca,'ydir','reverse') ; hold on ; title('f');
%   subplot(1,2,2); imagesc(x,x,real(cp)) ; axis image; colormap jet(256) ; colorbar ; set(gca,'ydir','reverse') ; hold on ; title('cp');
  
  clear cond ellipse1_* ellipse2* ev ew ex ez halfc_* k ilmax ixmin izmax izmin eh

end
