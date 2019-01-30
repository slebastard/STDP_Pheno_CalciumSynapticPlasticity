function S = poissonSampl(ns, lambda)
%POISSONSAMPL Summary of this function goes here
%   Detailed explanation goes here
  S = zeros(ns,1);
  for i = 1:ns    
      k = 0;
      produ = rand;
      while produ >= exp(-lambda)
          produ = produ*rand;
          k = k+1;
      end
      S(i) = k;
  end
end

