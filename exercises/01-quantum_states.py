"""
The purpose of this exercise is to familiarize you with some common ways to work
with our particular formulation of quantum states. Here, we'll go over how to
construct quantum states.
"""

#== Part 1: Making states ======================================================
print("\nPART 1")
# In this Python interface, the behaviors of quantum states will resemble (as
# closely as possible) those of the non-Python quantum states that we've already
# been working with. Let's start by making a set of basis states. All of this
# quantum business is not part of Python's standard loadout, so first we have to
# tell Python where to find it. We do this with an `import` statement, shown
# below:

from lib import *

# The exact details of what's going on here aren't so important, but you should
# know that you need this in any program you write here that does anything with
# quantum states.
#
# Now to make the actual basis. The general pattern of how this is accomplished
# is by starting with a `BasisState` and converting it to a regular `StateVec`,
# which looks like this:

my_basis_state = BasisState("label").to_statevec()

# Just one basis state is pretty boring, though, so we'll make some more:

a = BasisState("a").to_statevec()
b = BasisState("b").to_statevec()
c = BasisState("c").to_statevec()

# and just like that, we have three basis states |a>, |b>, and |c>! Try
# printing them out:

# ----- write your code here -----



# --------------------------------

# You may notice that every one of these are kets, shown by how they print out
# as `|...>` instead of `<...|`. So then how do we get bra versions of these
# states? The way to do this is by using the `.hc()` method. For example, if you
# wanted the bra versions of all of the states above, you'd write them as
#
# a.hc()
# b.hc()
# c.hc()
#
# Try printing them out:

# ----- write your code here -----





# --------------------------------



#== Part 2: Doing math with states =============================================
print("\nPART 2")
# Now let's try doing some math with the states from Part 1. First, let's
# redefine them:

a = BasisState("a").to_statevec()
b = BasisState("b").to_statevec()
c = BasisState("c").to_statevec()

# First, we importantly assume that every basis state created this way is
# orthonormal to each other. That is, we should expect that <a|a> == 1,
# <a|b> == 0, and <a|c> == 0. We can compute these products like so:

a.hc() * a
a.hc() * b
a.hc() * c

# Try wrapping these products above in `print` statements, and check to see that
# you get what you expect.
#
# But that's only a third of the story, though, since we have three basis states
# that we'd like to check in the same way. We've already checked all the
# relations for |a>, so now let's check those for |b> and |c>: compute and print
# out the products
#
# <b|a>
# <b|b>
# <b|c>
# <c|a>
# <c|b>
# <c|c>

# ----- write your code here -----







# --------------------------------

# We can also add states together. Let's say we wanted to make an arbitrary
# state
#
# |my_state> = 3 |a> + (1+i) |b> - 4i |c>
#
# We can do this with the following (remember that i = 1j in python):

my_state = 3*a + (1+1j)*b - 4j*c
#print(my_state)

# Try un-commenting the print statement above to make sure that it's what we
# expect it to be.
#
# Now for some interesting bits: what happens when we try multiplying `my_state`
# with some bras? Compute and print out the following products:
#
# <a|my_state>
# <b|my_state>
# <c|my_state>

# ----- write your code here -----




# --------------------------------

# We'll explore this more in the next exercise.

