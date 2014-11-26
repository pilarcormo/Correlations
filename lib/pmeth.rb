#encoding: utf-8

require 'prime'

class Integer
	#returns all positive factors of an Integer
	def factors
    	(1..self).select { |n| (self % n).zero? } 
	end
  
	#returns true if self is a prime number, false otherwise
	def prime?
		Prime.prime?(self)
	end
end


class PMeth

	# Returns a random integer that array can be divided by to get another integer (other than the array length itself)
	# by getting all the factors of the array length and ditching the array length
	def self.division(array)
		array.length.factors.select {|n| n != array.length}.sample
	end

	# Returns a new permutation that has had a randomly sized sub-section re-ordered by shuffle
	def self.chunk_mutate(permutation)
		if permutation.length.prime? # if there are a prime number of objects in the permutation
			ig = rand(permutation.length) # choose a random object to ignore - to add back at its original index after mutation
			ig_obj = permutation[ig] # save the object
			permutation.delete_at(ig)
		end
		x = 0 # x is the randomly chosen size of chunks that permutation will be split into
		until x >= 2 # the chunk to be re-ordered must have at least 2 objects
			x = division(permutation)
		end
		sliced = permutation.each_slice(x).to_a # permutation is sliced into chunks of size x
		e = rand(sliced.length) # one of the chunks is chosen at random...
		chunk = sliced[e]
		until sliced[e] != chunk
			sliced[e] = sliced[e].shuffle # ... and the objects within are shuffled
		end
		mutant = sliced.flatten
		if ig != nil
			mutant.insert(ig, ig_obj)
		end
		return mutant
	end

	# Returns a new permutation where two of the objects, chosen at random, have swapped positions (indices)
	def self.swap_mutate(permutation)
		a = b = x = y = 0
		until a != b
			x = rand(permutation.length) # randomly choose two indices x and y...
			y = rand(permutation.length)
			a = permutation[x] # ... and call the objects at these indices a and b
			b = permutation[y]
		end
		mutant = permutation.dup # create a new permutation...
		mutant[x] = b # ... with object b at index x...
		mutant[y] = a # ... and object a at index y
		return mutant
	end

	# Returns a new permutation where an object, chosen at random, has swapped position with one of it's immediate neighbours
	def self.adjacent_swap(permutation)
		x = rand(permutation.length) # randomly choose an index
		coin = rand(2) # a 50% chance for an object to swap with it's left or right neighbour
		if x == 0 # object at index 0 of array can only swap with index 1
			y = x + 1
		elsif x == permutation.length-1 # object at index -1 of array can only swap with index -2
			y = x - 1
		elsif coin == 0
			y = x + 1
		elsif coin == 1
			y = x - 1
		end
		mutant = permutation.dup
		mutant[x] = permutation[y]
		mutant[y] = permutation[x]
		return mutant
	end

	# Returns a permutation whose objects are ordered partly like parent_1 permutation, and partly like parent_2 permutation
	def self.recombine(parent_1, parent_2)
		x = division(parent_1) # the randomly chosen size of chunks that permutations will be split into
		if parent_1.length.prime? # to compensate for permutations with a prime number of objects:
			ig = rand(parent_1.length) # choose a random object to ignore - to add back at its original index after mutation
			p2_ig = parent_2.index(parent_1[ig]) # choose the same object from parent_2
			parent_1_reduced, parent_2_reduced = parent_1.dup, parent_2.dup # then create duplicates of the parent arrays...
			parent_1_reduced.delete_at(ig); parent_2_reduced.delete_at(p2_ig) # .. and remove the ignored object from these duplicates
			x = division(parent_1_reduced) # choose a new chunk size for reduced parent permutations, that no longer have a prime number of objects
			p1s, p2s = parent_1_reduced.each_slice(x).to_a, parent_2_reduced.each_slice(x).to_a # slice the reduced parent permutations into chunks of size x
		else
			p1s, p2s = parent_1.each_slice(x).to_a, parent_2.each_slice(x).to_a # if permutation lengths are non-prime, just slice the parent permutations into chunks of size x
		end
		chosen = rand(p2s.length) # choose a chunk to have parent_2 ordered objects in child permutation
		child = p1s.flatten.dup # un-modified child permutation to accept chunk from parent_2 (and possibly ignored object)
		p1s[chosen].each do |i| # place each object in chosen parent_1 chunk into the index it's corresponding object (from parent_2) occupies in parent_1:
			index = p1s[chosen].index(i) # the index of each object in the chosen parent_1 chunk...
			p2_object = p2s[chosen][index] # ... the object at that index in the chosen parent_2 chunk...
			p1_index = p1s.flatten.index(p2_object) # ... the index of that object (from parent_2 chunk) in parent_1
			unless p2s[chosen].include?(p1s[chosen][index]) # unless the chosen chunk from parent_2 includes objects from the chosen parent_1 chunk...
				child[p1_index] = p1s[chosen][index] # ... the object from the chosen chunk in parent_1 is added to p1_index in child
				chosen_index = p1s.flatten.index(p1s[chosen][index]) # index of object from parent_1's chosen chunk in parent_1
				child[chosen_index] = p2s[chosen][index] # swapping the indices of objects in chunks from parents, to give their indices in child
			end
		end
		if ig != nil # if the permutations are a prime number in length, we now add the ignored object from earlier into child
			if p2s[chosen].include?(parent_2[p2_ig]) # add the ignored object from parent_2 if it's in the chosen chunk...
				child.insert(p2_ig, parent_2[p2_ig])
			else
				child.insert(ig, parent_1[ig]) # ...otherwise add the ignored object from parent_1
			end
		end
		return child
	end
end