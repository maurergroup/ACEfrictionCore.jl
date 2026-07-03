


import Base:   ==
import ACEfrictionCore: read_dict, write_dict, 
       transform, transform_d, transform_dd, inv_transform
       
export polytransform, morsetransform, agnesitransform

abstract type DistanceTransform end        

# ----- new transforms implementation 
import ACEfrictionCore: λ 

polytransform(p, r0) = λ("r -> ((1+$r0)/(1+r))^$p")
@deprecate PolyTransform(p, r0) polytransform(p, r0)

idtransform() = λ("r -> r")
@deprecate IdTransform() idtransform()

morsetransform(lambda, r0) = λ("r -> exp(- $lambda * (r / $r0 - 1))")
@deprecate MorseTransform(lambda, r0) morsetransform(lambda, r0)

agnesitransform(r0, p=2, a=(p-1)/(p+1)) = λ("r -> 1 / (1 + $a * (r / $r0)^$p)")
@deprecate AgnesiTransform(args...) agnesitransform(args...)


# ------------------------------------------------------
#  implementation of inverse transform 

import Roots: find_zero
function inv_transform(trans, x)
   # solve trans(r) = x
   r = find_zero(r -> trans(r) - x, 1.0) 
   if !(trans(r) ≈ x)
      @warn("inv_transform via find_zero didn't find a good solution")
   end
   return r 
end 


# ------------------------------------------------------
# generic ad codes for distance transforms 

import ACEfrictionCore:  evaluate

evaluate(trans::DistanceTransform, r::Number) = transform(trans, r)