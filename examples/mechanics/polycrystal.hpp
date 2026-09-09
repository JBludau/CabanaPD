#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>

#include <Kokkos_Core.hpp>

// Get flat index into N-dim array
template <std::size_t n>
int indexND( const std::array<int, n>& index, const std::array<int, n>& shape )
{
    int outIndex = 0;
    int stride = 1;

    for ( int axis = n - 1; axis >= 0; --axis )
    {
        outIndex += index[axis] * stride;
        stride *= shape[axis];
    }

    return outIndex;
}

// Check if N-dim array multi-index is valid
template <std::size_t n>
bool isValid( const std::array<int, n>& idx, const std::array<int, n>& shape )
{
    for ( std::size_t i = 0; i < n; ++i )
    {
        if ( idx[i] < 0 || idx[i] >= shape[i] )
        {
            return false;
        }
    }
    return true;
}

// Recursive helper function for makeNeighborRelativeIndices
template <std::size_t n>
void makeRelativeIndicesRecursive(
    std::size_t axis, std::array<int, n>& curIndex,
    std::vector<std::array<int, n>>& outRelativeIndices )
{
    int maxDist = std::ceil( std::sqrt( static_cast<double>( n ) ) );
    for ( int i = -maxDist; i <= maxDist; ++i )
    {
        curIndex[axis] = i;
        if ( axis < n - 1 )
        {
            makeRelativeIndicesRecursive( axis + 1, curIndex,
                                          outRelativeIndices );
        }
        else
        {
            outRelativeIndices.push_back( curIndex );
        }
    }
}

// Get N-dim array neighbor relative multi-indices
template <std::size_t n>
void makeNeighborRelativeIndices(
    std::vector<std::array<int, n>>& outRelativeIndices )
{
    std::array<int, n> curIndex = {};
    makeRelativeIndicesRecursive( 0, curIndex, outRelativeIndices );
}

template <std::size_t n>
double distSquared( const std::array<double, n>& a,
                    const std::array<double, n>& b )
{
    double r2 = 0.0;
    for ( std::size_t axis = 0; axis < n; ++axis )
    {
        r2 += ( a[axis] - b[axis] ) * ( a[axis] - b[axis] );
    }
    return r2;
}

// n-dimensional Poisson disc sampling in a rectangular prism domain
template <std::size_t n, class RNGType>
void poissonDiscSampling( const std::array<double, n>& extent, double r, int k,
                          std::vector<std::array<double, n>>& outPoints,
                          std::array<int, n>& outGridShape, double& outCellSize,
                          RNGType& gen )
{
    // Precompute reused values for sampling
    double rn = std::pow( r, static_cast<double>( n ) );
    double rn2 = std::pow( 2.0 * r, static_cast<double>( n ) );

    // Calculate shape of grid
    double cellSize = r / std::sqrt( static_cast<double>( n ) );
    std::array<int, n> gridShape = {};
    int totalCells = 1;
    for ( std::size_t axis = 0; axis < n; ++axis )
    {
        double axisCells = std::ceil( extent[axis] / cellSize );

        // Add an extra cell in each direction just in case
        gridShape[axis] = 1 + static_cast<int>( axisCells );
        totalCells *= gridShape[axis];
    }

    // Initialize empty grid (as flat array)
    std::vector<int> grid( totalCells, -1 );

    // Calculate grid search relative indices
    std::vector<std::array<int, n>> nbrIndicesRel;
    makeNeighborRelativeIndices( nbrIndicesRel );

    // Choose first point
    std::uniform_real_distribution<double> coordDist( 0.0, 1.0 );
    auto coordGen = std::bind( coordDist, gen );

    std::array<double, n> x0 = {};
    std::array<int, n> idx0 = {};
    for ( std::size_t axis = 0; axis < n; ++axis )
    {
        x0[axis] = coordGen() * extent[axis];
        idx0[axis] = static_cast<int>( std::floor( x0[axis] / cellSize ) );
    };

    grid[indexND( idx0, gridShape )] = 0; // Store x0 location in grid
    outPoints.push_back( x0 );

    // Initialize active set
    std::vector<int> activeSet = { 0 };

    while ( activeSet.size() > 0 )
    {
        std::uniform_int_distribution<int> pointDist( 0, activeSet.size() - 1 );
        int seedIndexInActive = pointDist( gen );
        int seedIndex = activeSet[seedIndexInActive];
        std::array<double, n> seedX = outPoints[seedIndex];

        bool addedPoint = false;

        for ( int i = 0; i < k; ++i )
        {
            std::array<double, n> x = {};
            std::array<int, n> idx = {};

            // Use inverse method to sample distance from seed
            double xR = std::pow( rn + coordGen() * ( rn2 - rn ),
                                  1.0 / static_cast<double>( n ) );

            // Sample direction cosine angles uniformly in [0, pi] to get
            // direction
            bool failed = false;
            for ( std::size_t axis = 0; axis < n; ++axis )
            {
                x[axis] =
                    seedX[axis] + xR * std::cos( coordGen() * CabanaPD::pi );
                idx[axis] =
                    static_cast<int>( std::floor( x[axis] / cellSize ) );
                if ( x[axis] < 0 || x[axis] >= extent[axis] )
                {
                    failed = true;
                    break;
                }
            }
            if ( failed )
            {
                continue;
            }

            // Check if point is too close to any existing points
            std::array<int, n> checkIdx = idx;
            bool isClose = false;

            for ( const std::array<int, n>& relIdx : nbrIndicesRel )
            {
                for ( std::size_t axis = 0; axis < n; ++axis )
                {
                    checkIdx[axis] = idx[axis] + relIdx[axis];
                }

                if ( !isValid( checkIdx, gridShape ) )
                {
                    continue;
                }

                int checkPoint = grid[indexND( checkIdx, gridShape )];
                if ( checkPoint == -1 )
                {
                    continue;
                }

                const std::array<double, n>& checkX = outPoints[checkPoint];

                if ( distSquared( checkX, x ) > r * r )
                {
                    continue;
                }

                isClose = true;
                break;
            }

            // New point is not too close to any existing, so add it
            if ( !isClose )
            {
                outPoints.push_back( x );
                activeSet.push_back( outPoints.size() - 1 );
                grid[indexND( idx, gridShape )] = outPoints.size() - 1;
                addedPoint = true;
                break;
            }
        }
        // If no point could be generated farther than r from existing points,
        // remove seed point from active set
        if ( !addedPoint )
        {
            activeSet.erase( activeSet.begin() + seedIndexInActive );
        }
    }

    outGridShape = gridShape;
    outCellSize = cellSize;
}

// Generates numGrains random polycrystal grain centers using
// Poisson disc sampling
void getPolycrystalGrains( const std::array<double, 3>& extent,
                           double grainSize,
                           std::vector<std::array<double, 3>>& outLocations,
                           std::array<int, 3>& outGridShape,
                           double& outCellSize )
{
    // Initialize RNG
    std::minstd_rand baseRng;
    baseRng.seed( 12345 );
    std::seed_seq randomSeed{ baseRng(), baseRng(), baseRng(), baseRng(),
                              baseRng(), baseRng(), baseRng(), baseRng() };
    std::mt19937 gen( randomSeed );

    // Generate random, evenly-spaced candidate grain locations
    poissonDiscSampling( extent, grainSize, 30, outLocations, outGridShape,
                         outCellSize, gen );
}

template <typename GrainGridType, typename GrainPosType>
struct FindClosestGrainFunctor
{
    GrainGridType grainGrid;
    GrainPosType grainPos;
    double minX;
    double minY;
    double minZ;
    int sizeX;
    int sizeY;
    int sizeZ;
    double grainGridCellSize;

    FindClosestGrainFunctor( GrainGridType grainGrid, GrainPosType grainPos,
                             const std::array<double, 3>& low_corner,
                             const std::array<int, 3>& grainGridShape,
                             double grainGridCellSize )
        : grainGrid( grainGrid )
        , grainPos( grainPos )
        , minX( low_corner[0] )
        , minY( low_corner[1] )
        , minZ( low_corner[2] )
        , sizeX( grainGridShape[0] )
        , sizeY( grainGridShape[1] )
        , sizeZ( grainGridShape[2] )
        , grainGridCellSize( grainGridCellSize )
    {
    }

    KOKKOS_FUNCTION int operator()( double x, double y, double z ) const
    {
        int gx = static_cast<int>(
            Kokkos::floor( ( x - minX ) / grainGridCellSize ) );
        int gy = static_cast<int>(
            Kokkos::floor( ( y - minY ) / grainGridCellSize ) );
        int gz = static_cast<int>(
            Kokkos::floor( ( z - minZ ) / grainGridCellSize ) );

        // Search for any grain in shells of increasing size
        int shellRadius = 0;
        int guessIndex = -1;
        while ( guessIndex == -1 )
        {
            // X planes
            for ( int sy = Kokkos::max( gy - shellRadius, 0 );
                  sy <= Kokkos::min( gy + shellRadius, sizeY - 1 ); ++sy )
            {
                for ( int sz = Kokkos::max( gz - shellRadius, 0 );
                      sz <= Kokkos::min( gz + shellRadius, sizeZ - 1 ); ++sz )
                {
                    int sx = Kokkos::max( gx - shellRadius, 0 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sy = Kokkos::min( gy + shellRadius, sizeY - 1 ) + 1;
                        break;
                    }

                    sx = Kokkos::min( gx + shellRadius, sizeX - 1 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sy = Kokkos::min( gy + shellRadius, sizeY - 1 ) + 1;
                        break;
                    }
                }
            }

            if ( guessIndex > 0 )
                break;

            // Y planes
            for ( int sx = Kokkos::max( gx - shellRadius, 0 );
                  sx <= Kokkos::min( gx + shellRadius, sizeX - 1 ); ++sx )
            {
                for ( int sz = Kokkos::max( gz - shellRadius, 0 );
                      sz <= Kokkos::min( gz + shellRadius, sizeZ - 1 ); ++sz )
                {
                    int sy = Kokkos::max( gy - shellRadius, 0 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sx = Kokkos::min( gx + shellRadius, sizeX - 1 ) + 1;
                        break;
                    }

                    sy = Kokkos::min( gy + shellRadius, sizeY - 1 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sx = Kokkos::min( gx + shellRadius, sizeX - 1 ) + 1;
                        break;
                    }
                }
            }

            if ( guessIndex > 0 )
                break;

            // Z planes
            for ( int sy = Kokkos::max( gy - shellRadius, 0 );
                  sy <= Kokkos::min( gy + shellRadius, sizeY - 1 ); ++sy )
            {
                for ( int sx = Kokkos::max( gx - shellRadius, 0 );
                      sx <= Kokkos::min( gx + shellRadius, sizeX - 1 ); ++sx )
                {
                    int sz = Kokkos::max( gz - shellRadius, 0 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sy = Kokkos::min( gy + shellRadius, sizeY - 1 ) + 1;
                        break;
                    }

                    sz = Kokkos::min( gz + shellRadius, sizeZ - 1 );
                    guessIndex = grainGrid( sx, sy, sz );
                    if ( guessIndex >= 0 )
                    {
                        sy = Kokkos::min( gy + shellRadius, sizeY - 1 ) + 1;
                        break;
                    }
                }
            }

            ++shellRadius;
        }

        // Now find the nearest grain center
        double guessDX = grainPos( guessIndex, 0 ) - x;
        double guessDY = grainPos( guessIndex, 1 ) - y;
        double guessDZ = grainPos( guessIndex, 2 ) - z;
        double dist = Kokkos::sqrt( guessDX * guessDX + guessDY * guessDY +
                                    guessDZ * guessDZ );
        int gridSearchDist =
            static_cast<int>( Kokkos::ceil( dist / grainGridCellSize ) );

        double closestDist = dist;
        int closestIndex = guessIndex;
        for ( int sx = Kokkos::max( gx - gridSearchDist, 0 );
              sx <= Kokkos::min( gx + gridSearchDist, sizeX - 1 ); ++sx )
        {
            for ( int sy = Kokkos::max( gy - gridSearchDist, 0 );
                  sy <= Kokkos::min( gy + gridSearchDist, sizeY - 1 ); ++sy )
            {
                for ( int sz = Kokkos::max( gz - gridSearchDist, 0 );
                      sz <= Kokkos::min( gz + gridSearchDist, sizeZ - 1 );
                      ++sz )
                {
                    int testIndex = grainGrid( sx, sy, sz );
                    double testDX = grainPos( testIndex, 0 ) - x;
                    double testDY = grainPos( testIndex, 1 ) - y;
                    double testDZ = grainPos( testIndex, 2 ) - z;
                    double testDist = Kokkos::sqrt(
                        testDX * testDX + testDY * testDY + testDZ * testDZ );
                    if ( testDist < closestDist )
                    {
                        closestDist = testDist;
                        closestIndex = testIndex;
                    }
                }
            }
        }

        return closestIndex;
    }
};
