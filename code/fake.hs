-- Fake.hs

-- Creates fake data sets for our 3D scanner


{- 

We have two actors: the scanner and the object being scanned. The
scanner should have no information about the object, the object should
only reflect a measurement based on the geometry the scanner is trying
to reach.

Space is defined such that all useful z-positions are positive, the
bottom of the object should touch the plane z = 0. Theta positions
indicate the line of sight of the scanning head and have 2*pi
symmetry.

-}

import Numeric.GSL.Minimization

type Window = Double
type ZPos  = Double 
type ThPos = Double 
type ScanPos = (Window, ZPos, ThPos)
type ScanRange = (ZPos, -- zmin
                  ZPos, -- zmax
                  Double, -- zrez
                  Double) -- threz
type SpaceRange = (ZPos, -- zmin
                   ZPos, -- zmax
                   Double, -- zrez
                   Double) -- xrez

class Seeable a where
    measure :: a -> ScanPos -> Double

-- |Converts SpaceRange coordinates to ScanRange coordinates by
-- converting a x resolution to a necessary theta resolution, based on
-- a worst case scenario for the scanned object. 
scanCoord :: Double -- largest major elliptical axis
          -> Double -- largest minor elliptical axis
          -> SpaceRange
          -> ScanRange
scanCoord a b (z1, z2, zr, xr) = let (thr:_, _) = minfn (\th -> xr - error th)
                                 in (z1, z2, zr, thr)
    where error th = 2 * abs (r 0 - r th)
          r th = a * b / sqrt ( (a * sin(th))^2 + (b * cos(th))^2 )
          minfn f = 
              minimize NMSimplex2 1e-3 100 [1] (\[x] -> f x) [0]
scanCoord0 :: SpaceRange -> ScanRange
scanCoord0 = scanCoord 4 (4*1.15)
    
-- |scans with points uniformly distributed across (z,th) space. It's
-- not a realistic scan pattern, but it is a solid first test.
fakeScan :: Window -> ScanRange -> (ScanPos -> Double) -> [(ScanPos, Double)]
fakeScan w (zmin, zmax, zrez, threz) f = zip samples (map f samples)
    where samples = [(w, z, th) | z <- zs, th <- ths]
              where zs  = enumFromThenTo zmin (zmin + zrez) zmax
                    ths = enumFromThenTo 0 threz (2*pi)


data Sphere = Sphere Double -- Radius

instance Seeable Sphere where
    measure (Sphere r) (w, z, th) | 2*r*z - z^2 >= 0 = w/2 - sqrt(2*r*z - z^2)

fakesphere = let range = (scanCoord0 (0, 6, 0.002, 0.002))
                 obj   = measure (Sphere 4)
             in fakeScan 6 range obj
