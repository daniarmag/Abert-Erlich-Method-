import Data.Complex
import Text.Printf
import Data.Time.Clock
import Debug.Trace (trace)

--A function that receives Num type arguements -> a coeff list
--Returns the deriavtive of coeff list 
derivative :: Num a => [a] -> [a]
derivative [] = []
derivative [p] = []
--The functional part of the function when the array length is > 1
derivative (p:ps) = p' : derivative ps
  where
    --For each element, takes p and multplies with length of the remaining ps list (equals to pwoer of p)
    p' = fromIntegral (length ps) * p

-- A function that creates a list of powers and a list of complex representations of coeff 
-- uses them to calculate p(x) = an*x^n+.....+a1*x+a0 (SLIDE 4).
polyval :: [Double] -> Complex Double -> Complex Double
polyval p x = sum (zipWith (*) (reverse (map toComplex p)) (powers x))
  where
    --when calling d - reversing the list to match the power list order.
    toComplex d = d :+ 0
    --Creates an infinite list of coeff powers. Uses lazy eval structure so 
    --only the necessary elements are computed in each function call.
    powers x' = iterate (*x') (1 :+ 0)

-- This function simply divides value of poly at x and poly' at x.
-- Two cases here to avoid overflow for very large numbers.
-- Receives two double arrays (p and p') and a complex double (x)
divide :: [Double] -> [Double] -> Complex Double -> Complex Double
divide p q x
  --SLIDE 8, different calculation to avoid the overflow problem
  | magnitude x <= 1 = polyval p x / polyval q x
  --here we can see the reverses part -> p(zk)/q(zk)=zk*(reversed_p(1/zk)/reversed_q(1_zk))
  | otherwise = x * (polyval (reverse p) (1 / x) / polyval (reverse q) (1 / x))

--Helper function for init roots, receives R and theta as double type and returns the complex representation
--(SLIDE 4)
eulerEquation :: Double -> Double -> Complex Double
eulerEquation r theta =
  let realPart = r * cos theta
      imagPart = r * sin theta
  in realPart :+ imagPart

--Calculates the initial roots. Entire implementation of equation from SLIDE 4 
initRoots :: [Double] -> [Complex Double]
initRoots p = [ eulerEquation r  (2 * pi * fromIntegral k / fromIntegral n) | k <- [0..n-2] ]
  where
    --normalizing the radius by diving the calculation with the first coeff
    --adding a safety epsilon to avoid division by zero
    r = 1 + maximum (map abs (drop 1 p)) / (abs (head p) + 1e-10)
    n = length p

--SLIDE 5, the sigma part.
sums :: Fractional a => [a] -> [a]
--Nested loop. Outer loop k, inside loop j.
--for each k, calculate the sum of 1/zk-zj where k!=j
sums s = [sum [1 / (zk - zj) | (zj, j) <- zip s [1..], j /= k] | (zk, k) <- zip s [1..]]

--SLIDE 5, entire equation.
--Receives coeff list, derivative list, roots list and returns offsets list (type complex double) 
getOffsets :: [Double] -> [Double] -> [Complex Double] -> [Complex Double]
getOffsets p pTag roots =
  --for each root in roots, apply divide function with p and ptag lists
  let numerator = map (divide p pTag) roots
      calcSigma = sums roots
      --for each n in numerator and s in calcSigma
      denominator = zipWith (\n s -> 1 - (n * s)) numerator calcSigma
      --for each element in the corresponding lists, divide the elemnts and put all in a list
      wList = zipWith (/) numerator denominator
  in wList

-- Combining all the functions to eventually implement the alg and return the final roots
abertErlich :: [Double] -> [Double] -> [Complex Double] -> Double -> Int -> [Complex Double]
--Private function for each iteration of the algorithm
abertErlich p pTag roots epsilon maxTries = abertErlich' roots 0
  where
    abertErlich' currentRoots tries =
      --w is a list of current offsets (full equation in SLIDE 5)
      let w = getOffsets p pTag currentRoots
          --substract each root and its offset to get the updated roots.
          roots' = zipWith (-) currentRoots w
          --the alg stops when each offset is smaller than the defined epsilon
          --to avoid checking each offset, we simply take the max one and checking if it's smaller than epsilon
          maxW = maximum (map magnitude w)
          debugMsg = trace ("maxW: " ++ show maxW) maxW
      in if debugMsg < epsilon || tries >= maxTries
          -- if the condition is met, return the final roots
           then roots'
           --recursive call of inner function abertErlic' as long as we did not reach convergence.
           else abertErlich' roots' (tries + 1)

--Helper function that returns the context of the file input as a list of doubles
readFileToList :: FilePath -> IO [Double]
readFileToList filePath = do
    -- read entire file to contents string
    contents <- readFile filePath
    -- divide contents into strings, cast each string to double and map them all to a list
    let numbers = map read (lines contents) :: [Double]
    return numbers

--This main function constructs the general flow of Abert Erlich algorithm 
--Also manages parameters that help control when the algorithm reaches an end
main :: IO ()
main = do
    start <- getCurrentTime
    let maxTries = 850
    let epsilon = 1e-5
    p <- readFileToList "poly_coeff(997).txt"
    let dp = derivative p
    let roots = initRoots p
    let abertErlichRoots = abertErlich p dp roots epsilon maxTries
    putStrLn "Found roots:"
    mapM_ (\root -> printf "%.3f%+.3fj\n" (realPart root) (imagPart root)) abertErlichRoots
    end <- getCurrentTime
    let duration = diffUTCTime end start
    putStrLn ( "The operation took: " ++ show duration ++ " seconds" )
