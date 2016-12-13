import numpy as np

from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.helper import findSetBits, sortPermutations

from work.probabilistic_transformations_for_inference.sampling import TensorQuadratureSampleGenerationStrategy, \
    ApproximateFeketeSampleGeneratorStrategy, LejaSampleGeneratorStrategy
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes

class TestEnvironmentSG(object):

    def buildSetting(self,
                     f,
                     params,
                     level,
                     gridType,
                     deg=1,
                     maxGridSize=1000,
                     isFull=False,
                     epsilon=1e-15,
                     adaptive=None,
                     adaptPoints=3,
                     uqSetting=None,
                     uqSettingRef=None,
                     knowledgeFilename=None):
        builder = ASGCUQManagerBuilder()

        builder.withParameters(params)\
               .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                      KnowledgeTypes.SQUARED,
                                      KnowledgeTypes.EXPECTATIONVALUE])\
               .useInterpolation()

        builder.defineUQSetting().withSimulation(f)

        if uqSetting is not None:
            build.useUQSetting(uqSetting)

        if uqSettingRef is not None and len(uqSettingRef) > 0:
            builder.withTestSet(uqSettingRef)\
                   .learnWithTest()

        if knowledgeFilename is not None:
            builder.withKnowledge(knowledgeFilename)

        samplerSpec = builder.defineSampler()
        gridSpec = samplerSpec.withGrid()
        gridSpec.withLevel(level).hasType(gridType)
        if deg > 1:
            gridSpec.withDegree(deg)
        if isFull:
            gridSpec.isFull()

        if adaptive is not None:
            # specify the refinement
            samplerSpec.withRefinement()\
                       .withAdaptThreshold(epsilon)\
                       .withAdaptPoints(adaptPoints)\
                       .withBalancing()

            refinement = samplerSpec.withRefinement().refineMostPromisingNodes()
            refinement.createAllChildrenOnRefinement()
            if adaptive == "simple":
                refinement.withSurplusRanking()
            elif adaptive == "exp":
                refinement.withExpectationValueOptimizationRanking()
            elif adaptive == "var":
                refinement.withVarianceOptimizationRanking()
            elif adaptive == "squared":
                refinement.withSquaredSurplusRanking()

            samplerSpec.withStopPolicy().withGridSizeLimit(maxGridSize)

        uqManager = builder.andGetResult()

        # update the stats, which are not stored with the knowledge
        if knowledgeFilename is not None:
            uqManager.recomputeStats()

        return uqManager

    def runSampler(self, uqManager, label, alabel, blabel):
        uqSetting = self.uqSettings[label]
        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, alabel, blabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # write the setting to file
        uqManager.uqSetting.writeToFile()


class ProbabilisticSpaceSGpp(object):
    
    def __init__(self, numDims):
        self.numDims = numDims
    
    def uniform(self, a=0, b=1):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.numDims):
            up.new().isCalled("x%i" % idim).withUniformDistribution(a, b)
        return up.andGetResult()

    def normal(self, mu=0, sigma=1, alpha=0.001):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.numDims):
            up.new().isCalled("x%i" % idim).withNormalDistribution(mu, sigma, 0.001)
        return builder.andGetResult()

    def multivariate_normal(self, mu=0, cov=None, a=0, b=1):
            # set distributions of the input parameters
        mu = np.array([mu] * numDims)
        if cov is None:
            # use standard values
            diag = np.eye(numDims) * 0.005
            offdiag = np.abs(np.eye(numDims) - 1) * 0.001
            cov = diag + offdiag
        # estimate the density
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        names = ", ".join(["x%i" for i in xrange(numDims)])
        up.new().isCalled(names).withMultivariateNormalDistribution(mu, cov, 0, 1)
        return builder.andGetResult()


class PCEBuilderHeat(object):

    def __init__(self, numDims):
        self.numDims = numDims


    def define_expansion(self, pce, expansion, degree_1d):
        if expansion == "full_tensor":
            pce.define_full_tensor_expansion(degree_1d)
        elif expansion == "total_degree":
            pce.define_isotropic_expansion(degree_1d, 1.0)
        else:
            raise AttributeError("expansion '%s' is not supported" % expansion)

    def define_full_tensor_samples(self, sample_type, rv_trans, expansion):
        return TensorQuadratureSampleGenerationStrategy("uniform", rv_trans, expansion)
    
    def define_approximate_fekete_samples(self, samples, pce, rv_trans):
        return ApproximateFeketeSampleGeneratorStrategy(samples, pce, rv_trans)
    
    def define_approximate_leja_samples(self, samples, pce, rv_trans):
        return LejaSampleGeneratorStrategy(samples, pce, rv_trans)

    def eval_samples(self, samples, rv_trans, f):
        trans_samples = rv_trans.map_from_canonical_distributions(samples)
        values = np.ndarray(trans_samples.shape[1])
        for i, sample in enumerate(trans_samples.T):
            values[i] = f(sample)
        return values

    def getSortedSobolIndices(self, pce):
        sobol_indices = pce.sobol_indices()
        indices = [findSetBits(i + 1) for i in xrange(len(sobol_indices))]
        indices, ixs = sortPermutations(indices, index_return=True)
        sobol_indices_dict = {}
        for index, i in zip(indices, ixs):
            sobol_indices_dict[index] = sobol_indices[i]
        return sobol_indices_dict
