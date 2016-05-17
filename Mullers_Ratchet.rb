include Math
require 'test/unit'
include Test::Unit::Assertions

class Ratchet

	def initialize
		@big_N = 2000
		@genome_size = 1000
		@gens = 2000
		@mutation_rate = 1.2
		@mutations_to_track = 0
		@outFile = $stdout
    	@outFile = open(ARGV[0], 'a') if ARGV[0]
    	@prng = Random.new
		@selection_coefficient = 0.02
	   	@timer = 0
	end

	def iterate
		@iteration = true
		1.times{time_passes(@big_N, @genome_size, @gens, @mutation_rate, @mutations_to_track, @prng, @selection_coefficient)}
		@outFile.close if ARGV[0]
	end

	# the output of this program is an array indexed by 
	# number of mutations, such that if ratchet[i] = j, 
	# then the last organism with i or fewer mutations 
	# was lost at generation j.
	def time_passes(big_N, genome_size, gens, mutation_rate, mutations_to_track, prng, selection_coefficient)
		classes = Array.new(genome_size) {0}
		classes[0] = big_N
		ratchet = Array.new(mutations_to_track) {0}
		min = 0
		fitnesses = fitness_calc(genome_size, selection_coefficient)
		mut_probs = mut_probs_calc(genome_size, mutation_rate)
		(0...gens).each do |i|
			t = find_t(classes, fitnesses)
			break if (t < 1)
			classes = mutate(big_N, classes, fitnesses, genome_size, mut_probs, prng, t)
			new_min = classes.index{|m| !(m == 0)}
			puts "#{@timer}\ti #{i} \tnew_min #{new_min}" if @iteration 
			if (new_min and new_min > min)
				((min)..(new_min)).each {|j| ratchet[j] = i}
				min = new_min
			end
			break if (new_min and new_min > mutations_to_track)
		end
		@outFile.write "\n #{@ratchet} \n"
		@timer += 1
		ratchet
	end

	# The selection model implemented is that of multiplicative deleterious 
	# mutations.  Each mutation is independent and decreases fitness by a set
	# proportion, the selection coefficient.  Thus an entity with 0 
	# mutations has a fitness of 1; with 3 mutations, it has a fitness of
	# (1-s)^3.  This subroutine creates a precalculated array of fitnesses
	# for the possible numbers of mutations in the genome.
	def fitness_calc(genome_size, selection_coefficient)
		fitnesses = Array.new(genome_size) {0}
		running_fitness = 1
		(0...genome_size).each do |i|
			fitnesses[i] = running_fitness
			running_fitness *=(1-selection_coefficient)
		end
		fitnesses
	end

	# 1-(@mutation_rate/@genome_size)is the probability that a given gene will 
	# not mutate.  The probability of 0 mutations in a genome is
	# the probability that all genes will not mutate, 
	# p = (1-(@mutation_rate/@genome_size))**@genome_size
	# The limit as @genome_size-> inf. = exp(-@mutation_rate)
	# The limit as @genome_size-> inf. for j mutations is e^-U U^j over j!, 
	# from the def. of 'e'.  I use the limit because it doesn't invite the 
	# long float errors of the exact expression.
	def mut_probs_calc(genome_size, mutation_rate)
		mut_probs = Array.new(genome_size) {0}
		p = exp(-mutation_rate)
		mut_probs[0] = p
		(1...genome_size).each do |i|
			p*= mutation_rate/i.to_f # this is the probability of exactly j mutations per entity
			mut_probs[i] = p
		end
		mut_probs
	end

	# Each class represents the number of entities with k mutations.  
	# mutate iterates through the classes.  For each class k, it looks at every
	# class k-j >= 0, such that a member of class k-j could acquire j
	# mutations and enter class k.  The number of entities from k-j that are
	# transferred to k are computed using the number of entities in class k-j, 
	# their fitness, and the probability of exactly j mutations.

	def mutate(big_N, classes, fitnesses, genome_size, mut_probs, prng, t)
		new_classes = Array.new(genome_size){0}
		#break if t < 1
		classes.each_index do |k|
			p = 0	
			(1..(k+1)).each do |j|
				j -=1
 				p += classes[k-j] * fitnesses[k-j] * mut_probs[j] / t
			end
			p *= big_N
			new_classes[k] = poisson_calc(p, prng)  #This version uses only poisson, instead of the 
											  #normal approximation used in Haigh 1978
		end
		new_classes
	end

	def find_t(classes, fitnesses)
		t = 0
		classes.each_with_index do |value, i|
			t += (value*fitnesses[i])
		end
		t
	end

	# This is Knuth's algorithm for Poisson RNG.  
	# The sensible thing would be to do some kind of precalculation
	#  to reduce computation time
	def poisson_calc(expected, prng) #expected number of occurences
		l = Math.exp(-expected) #(P|X = 0)    
		k = 0.0 #(number of occurences)
		p = 1.0 #The product of the urandoms
		k = 1.0 unless p > l # need because 1.0 !> 1.0
		while p > l   
			u = prng.rand
			p *= u
			k += 1 
		end
	    (k - 1)
	end

end

class Test_Case < Test::Unit::TestCase
	def test_poisson_calc
		prng = Random.new(1234)
		input = 4
		expectedResult = 5
		actualResult = Ratchet.new.poisson_calc(input, prng)	
		assert_equal(expectedResult, actualResult)	
	end	

	def test_find_t
		genome_size = 20
		classes = [600, 300, 150, 75]
		selection_coefficient = 0.01
		fitnesses = Ratchet.new.fitness_calc(genome_size, selection_coefficient)
		expectedResult = 1116.787425
		actualResult = Ratchet.new.find_t(classes, fitnesses)
		assert_equal(expectedResult, actualResult)
	end

	def test_mutate
		genome_size = 20
		classes = [600, 300, 150, 75]
		selection_coefficient = 0.01
		mutation_rate = 1
		fitnesses = Ratchet.new.fitness_calc(genome_size, selection_coefficient)
		mut_probs = Ratchet.new.mut_probs_calc(genome_size, mutation_rate)
		t = Ratchet.new.find_t(classes, fitnesses)
		big_N = 1000
		prng = Random.new(1234)
		expectedResult = [216.0, 295.0, 211.0, 174.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
		actualResult = Ratchet.new.mutate(big_N, classes, fitnesses, genome_size, mut_probs, prng, t)
		assert_equal(expectedResult, actualResult)	
	end


	def test_mut_probs_calc
		mutation_rate = 1
		genome_size = 20
		actualResult = Ratchet.new.mut_probs_calc(genome_size, mutation_rate)
		expectedResult = [0.36787944117144233, 0.36787944117144233, 0.18393972058572117, 0.061313240195240384, 0.015328310048810096, 0.0030656620097620196, 0.0005109436682936699, 7.299195261338141e-05, 9.123994076672677e-06, 1.0137771196302974e-06, 1.0137771196302975e-07, 9.216155633002704e-09, 7.68012969416892e-10, 5.907792072437631e-11, 4.2198514803125934e-12, 2.8132343202083955e-13, 1.7582714501302472e-14, 1.0342773236060278e-15, 5.745985131144598e-17, 3.0242027006024198e-18]
		assert_equal(expectedResult, actualResult)
	end

	def test_fitness_calc
		genome_size = 20
		selection_coefficient = 0.01
		expectedResult = [1, 0.99, 0.9801, 0.9702989999999999, 0.96059601, 0.9509900498999999, 0.9414801494009999, 0.9320653479069899, 0.92274469442792, 0.9135172474836407, 0.9043820750088043, 0.8953382542587163, 0.8863848717161291, 0.8775210229989678, 0.8687458127689781, 0.8600583546412883, 0.8514577710948754, 0.8429431933839266, 0.8345137614500874, 0.8261686238355865]
		actualResult = Ratchet.new.fitness_calc(genome_size, selection_coefficient)
		assert_equal(expectedResult, actualResult)
	end

	def test_time_passes
		prng = Random.new(1234)
		genome_size = 20
		big_N = 200
		gens = 20
		mutations_to_track = 10
		selection_coefficient = 0.01
		mutation_rate = 1
    	expectedResult = [3, 5, 6, 10, 10, 11, 11, 13, 13, 16, 16, 16]
    	actualResult = Ratchet.new.time_passes(big_N, genome_size, gens, mutation_rate, mutations_to_track, prng, selection_coefficient)
    	assert_equal(expectedResult, actualResult)
	end

end

Ratchet.new.iterate # starts simulations

__END__

This program codes the algorithm described in Haigh 1978 for calculating the
rate of Muller's Ratchet in asexual populations.  In a nutshell, the Ratchet, 
first suggested in Muller ?, is a process by which asexual populations 
suffer higher mutational loads than sexual populations.  The mechanism is that, 
in any population, the "least-loaded class" -- the set of organisms with the 
fewest mutations -- may go extinct in a given generation by chance.  Sexual 
populations can regenerate such a least-loaded class through sexual
recombination.  Asexuals cannot; each time a least-loaded class is lost, the 
ratchet clicks forward, and the asexuals are weaker.  

The speed of the ratchet is dependent on the size of the population, the 
rate of mutation, and the selection coefficient of the mutations.  This version of
the program is intended for evolutionary biologists learning Ruby, so the 
task of iterating over many values of N, U, and s is left as an exercise for the
student.



# The following notes are useful for translating from Haigh 1978 to the Perl code.  

# p_k(t+1) = Sum from j = 0 to k of X_k-j(t) * (1-s)**(k-j) * [e**-U * U**j / j!] / T_1(t) 

# e**-U * U**j / j! = e**-U * the product from i = 0 to j of U/j.
#this is read out of a vector @mut_probs created by "generate_mut_probs"
# such that $mut_probs[j] == e**-U * U**j /j!

# (1-s)**h = Product from i = 0 to h of (1-s)
#this is read out of a vector @fitnesses created by "generate_fitnesses"
#such that $fitnesses[h] == (1-s)**h

# T_n(t) = Sum from i = 0 to inf. of X_i(t) * (1-s)**(i*r)
# T_1(t) = Sum from i = 0 to inf. of X_i(t) * (1-s)**i
#this is generated through the "find_T" function, that iterates over $classes[i]*$fitness[$i]
# remember T_1(t) = N * W_bar(t) (size of population * population mean fitness)

"A simple generator for random numbers taken from a Poisson distribution 
is obtained using this simple recipe: if x1, x2, ... is a sequence of random 
numbers with uniform distribution between zero and one, 
k is the first integer for which the product x1 · x2 · ... · xk+1 < e-λ" 
(untraced citation, attributed to Knuth)



