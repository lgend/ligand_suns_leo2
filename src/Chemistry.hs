{- Copyright 2013 Gabriel Gonzalez

   This file is part of the Suns Search Engine

   The Suns Search Engine is free software: you can redistribute it and/or
   modify it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or (at your
   option) any later version.

   The Suns Search Engine is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   the Suns Search Engine.  If not, see <http://www.gnu.org/licenses/>.
-}

-- | Substructure matching code

{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Chemistry
    ( -- * Type
      ParseS(..)
    , evalParseS

      -- * Substructure matching
    , pMotif2
    ) where
import Debug.Trace (traceShow)
import AtomName (AtomName)
import Control.Applicative (Applicative, Alternative)
import Control.Error (justZ)
import Control.Monad (MonadPlus, forM_, guard)
import Control.Monad.Trans.Class (lift)
import Control.Monad.Trans.State (
    StateT(StateT), evalStateT, get, put )
import Control.Monad.Trans.Reader (ReaderT, runReaderT)
import qualified Control.Monad.Trans.Reader as ReaderT (ask)
import Data.Array ((!))
import Data.Set (Set, insert, member)
import qualified Data.Traversable as T
import qualified Data.Vector.Generic as G
import qualified Data.Vector as V
import Data.Vector.Unboxed (Vector)
import qualified Data.Vector.Storable as VS
import Data.Set (empty)
import Structure (Structure(Structure), bonds, deleteBond)

-- | A structure parser
--newtype ParseS a = ParseS { unParseS :: StateT Structure [] a }
--    deriving (Functor, Applicative, Monad, Alternative, MonadPlus)

type AtomIndex = Int
type Match = VS.Vector AtomIndex
type MatchesSoFar = Set Match

newtype ParseS a = ParseS { unParseS :: ReaderT Structure (StateT MatchesSoFar []) a }
    deriving (Functor, Applicative, Monad, Alternative, MonadPlus)

{-| Run a structure parser on a 'Structure' graph, returning a list of
    alternative solutions
-}

{-
evalParseS :: ParseS a -> ReaderT Structure (StateT MatchesSoFar []) a -> [a]
evalParseS p = evalStateT (unParseS p)
-}
evalParseS :: ParseS a -> Structure -> [a]
evalParseS parser structure = evalStateT (runReaderT (unParseS parser) structure) empty

{-| Match a bond from the current graph and remove it from the graph

    'pBond' returns a pair of matched vertices in the same order as the argument
    'AtomName's.
-}

pBond :: AtomName -> AtomName -> ParseS (Int, Int)
pBond name1 name2 = ParseS $ do
    Structure gr as <- ReaderT.ask  -- ask is a `ReaderT` command
    -- Find all indices that match the first AtomName
    i1 <- lift $ lift $ VS.toList $ VS.findIndices (== name1) as
    -- Find all neighbors of the first atom that match the second AtomName
    i2 <- lift $ lift $ filter (\i -> as VS.! i == name2) (gr ! i1)
    -- Return the bond and the updated graph
    return (i1, i2)	

pMotif :: Structure -> ParseS (VS.Vector Int)
pMotif (Structure bs as)
    -- Use the State monad to keep track of matches
  = (`evalStateT` (V.replicate (VS.length as) Nothing)) $ do

        -- foreach (i1, i2) in (bonds graph):
        forM_ (bonds bs) $ \(i1, i2) -> do

            -- Match the bond
            (i1', i2') <- lift $ pBond (as VS.! i1) (as VS.! i2)

            -- The match must be consistent with other matches
            matches    <- get
            let consistent i1 i1' = case (matches V.! i1) of
                    Nothing   -> True
                    Just iOld -> iOld == i1'
            guard $ consistent i1 i1' && consistent i2 i2'

            -- Update the match list
            put $ matches V.// [(i1, Just i1'), (i2, Just i2')]

        -- Return the final list of matches
        matchesFinal <- get
        fmap G.convert $ justZ $ T.sequence matchesFinal

pMotif2 :: Structure -> ParseS (VS.Vector Int)
pMotif2 structure = ParseS $ do
   match <- unParseS (pMotif structure) 
   matchesSoFar <- lift get
   guard (not (match `member` matchesSoFar))
   lift $ put $! insert match matchesSoFar
   return match
