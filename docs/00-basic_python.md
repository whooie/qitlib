# The Basics of Python

Python is an interpreted, high-level programming language known for its
flexibility and wide applicability in addition to its accessibility to new
programmers. This is a short description of some key features of Python as they
relate to this project. It is hardly comprehensive, however; if you'd like to
learn more about this language and how to use it, I suggest taking a look at the
[official beginner's guide][1] or this [series of recorded lectures from an MIT
course][2].

## 0. Structure of a Python program
The fundamental goal of any programming language is to describe a well-ordered
list of explicit actions for a computer to complete, rather like writing down a
recipe. For example, say you wanted to instruct someone on how to make a peanut
butter-and-jelly sandwich. How would you go about this? After taking a second or
two to jot it down, you might start with something like this:
```
1. Put two slices of bread on a plate;
2. Spread peanut butter over one slice;
3. Spread jelly over the other slice;
4. Put the two slices together;
5. Done.
```

This recipe is fine for pretty much anyone; if anything, it might even be too
detailed. But what if you had to write the recipe for an alien? Even assuming
the alien knew what "bread", "slices", "peanut butter", "jelly", and "plates"
are, it would still be left with many questions, like:
- Where do I get the bread, plate, peanut butter, and jelly?
- How much of the peanut butter and jelly should I use?
- How do I spread them over the slices of bread?
- What sides of the slices of bread should I use for this?
- How should I put the slices together at the end?
- What should I do if the bread isn't pre-sliced or I can't open the jars of
  peanut butter and jelly?

The key point to remember is that computers (playing the role of the alien here)
are very dumb, and can only accomplish complicated tasks through explicit,
detailed instruction. When writing a program, you should always keep in mind
that the computer running your program knows nothing except for what you tell it
and a base set of rules that come with the language.

Additionally, some programming languages require the computer to look over a
complete set of instructions (that is, the entire program) before actually
deciding what to do, but Python does not. Whereas one of these other languages
may for some reason be able to see, like you or I might, that the order in which
steps (2) and (3) in the PB&J recipe doesn't matter, Python will not. Instead,
it forces the computer to look at each line of the program and execute it
one-by-one. This means that the general structure of a program written in Python
should look broadly like this one, which implements a well-known algorithm for
finding all the prime numbers less than some number `N`:
```python
# define basic names and terms for the program at the top of the program.
# by default, Python doesn't know what square roots are, so we have to define it
# here using someone else's code
from math import sqrt

# we want this to find all the prime numbers that are less than or equal to N.
# we'll do this by looking at each number less than or equal to N and checking
# to see if it's divisible by any other number in the list. Any number that is
# not divisible by any other number in this list is prime.
def sieve_of_eratosthenes(N):

    # start with a list of all integers in the range [2, N + 1)
    candidates = [n for n in range(2, N + 1)]

    # do the following for each of the candidate numbers
    for n in candidates:

        # we don't need to check any numbers greater than √N
        if n > sqrt(N):
            break

        # we also don't need to check any numbers we've already eliminated
        if n == -1:
            continue

        # set all of n's multiples equal to -1
        for j in range(n**2 + n - 2, N - 1, n):
            candidates[j] = -1

    # any remaining numbers not equal to -1 is prime!
    primes = [n for n in candidates if n != -1]
    return primes

# print all primes up to 100
print(sieve_of_eratosthenes(100))
```

Don't worry about understanding all of this right away. We'll break it down in
the following sections while we go through some of the basic building blocks of
the Python language. As a starting point, though, you may have noticed some
liberal use of the `#` character -- this signals a 'comment' in code, and tells
Python to ignore whatever comes after it. Using comments to jot down brief notes
is a very useful way to keep track of what you're doing and generally make your
programs more understandable to other people.

## 1. Variables and data types
Nearly all programming languages differentiate between types of data. Different
languages have different standards for how its types are dealt with (Python is
particularly flexible in this regard), but they all require the programmer to
think about, for instance, the difference between the number `12` and the letter
`'a'`. This is partly to help the computer manage all the bits and bytes that
make up a program while it's running, and partly because some operations simply
don't make sense for some types. (What would it mean to compute the sum `12 +
'a'`?)

### Numbers
The first type that we'll deal with is the **number**. In Python, "number"
actually comprises four data types that can all mutually interact through
standard arithmetic operations (`+`, `-`, `*` `/`).

#### `int`
`int`s are how Python represents integer numbers (that is, all the whole numbers
from minus infinity to positive infinity). `int`s can add to, subtract from, and
multiply each other, but there's a bit of fine print when it comes to division.
In general, the division of two integers will not itself be another integer
(consider the operation `5 / 2`), so Python actually has two kinds of division,
denoted by `//` and `/`. The first kind is what you should use when you want to
divide two `int`s and get another `int` out: it computes the full division and
discards any non-integer remainder (`5 // 2` equals `2`). The second kind is
what you should use if you want the full answer as floating-point number, which
brings us to...

#### `float`
`float`s represent floating-point numbers, or numbers that may have a
non-integer part. `float`s can interact with all other numbers in the way you'd
normally expect, but can be written in two different but functionally equivalent
ways. The first is the standard way, e.g. `250.734`, and the second is with
scientific notation, e.g. `2.50734e2` as shorthand for 2.50734 \* 10^2 or, in
Python syntax, `2.50734 * 10**2`.

#### `complex`
`complex`s are actually two `float`s bundled together to represent the real and
imaginary parts of a single complex number. To use a `complex`, you can write it
using `j` notation, e.g. `1.5 + 0.4j`, `-2.6 + 1e6j`. When two complex numbers
are added, subtracted, multiplied, or divided together, the imaginary components
obey the usual rules, i.e. `1j * 1j` is equal to `-1`.

#### `bool`
`bool`s (short for Boolean values) represent truth values for logical
expressions, and can be either `True` or `False`. When interacting with the
other number types, `bool`s take on the values `False = 0` and `True = 1`. You
can also work with `bool`s through *boolean expressions*, which are like regular
arithmetic ones, except they return either `True` or `False` instead. For
example, say we want to perform some checks on the number `5` to test its value.
To this purpose, we can use Python's `==`, `!=`, `<`, `<=`, `>`, and `>=`
operators which are equivalent to asking whether a value is equal, not equal,
less than, less than or equal to, greater than, and greater than or equal to
another value. When used, each of these symbols will return either `True` or
`False`:
```python
5 == 8 # False
5 != 8 # True
5 <  8 # True
5 <= 8 # True
5 >  8 # False
5 >= 8 # False
5 == 5 # True
5 != 5 # False
5 >= 5 # True
```

### Variables
What is a variable? Simply put, it's a way to store a value in the computer's
memory so that it can be referred back to at a later time. At its most basic, a
variable is just a name you give to some value, just like in regular math. In
Python, values are assigned to names with the `=` operator like so:
```python
my_num1 = 5
my_num2 = 8
```

Now we have the values `5` and `8` stored in the variable names `my_num1` and
`my_num2`, respectively. Notice that, unlike the usual your usual math homework,
our variable names can be more than one letter, and include special characters
like `_` and even numbers. In fact, anything not beginning with a number or
containing any punctuation other than `_` is a valid variable name in Python, so
you could even set `我 = 5` and `你 = 8` if you wanted (for the sake of
readability, though, please don't)!

Now that we've assigned values to these names, we can use them in any situation
where we might otherwise use the values they refer to. For example, we can
perform all the arithmetic operations on them, and even assign them to other
variables:
```python
my_num1 + my_num2
my_num1 * my_num3
my_num1 < my_num2
...
my_num3 = my_num1
my_num4 = my_num3 + my_num2
my_bool = my_num1 == my_num2
```

Notice that `=` (single equals) is *variable assignment*, while `==` (double
equals) is the *equality operator*. Being careful about which of these is used
and where will save you a few headaches down the line.

### Lists
In the example above, we created a lot of variables with similar names.
`my_num1`, `my_num2`, `my_num3`, and `my_num4` all refer to numbers that are
related to each other, so it would make sense to collect them all together so
that we don't have to keep track of four independent names in the rest of the
program, wouldn't you say? `list`s are the primary way that Python provides to
to this. Specifically, `list`s are ordered, untyped collections of data, all
grouped under one variable name, and are denoted with `[...]` square brackets.

### Strings

### Dictionaries

## 2. Logic and control flow

### `if` statements

### `for` loops

### `while` loops

## 3. Functions

## 4. General tips for writing Python

[1]: https://wiki.python.org/moin/BeginnersGuide
[2]: https://www.youtube.com/playlist?list=PLUl4u3cNGP63WbdFxL8giv4yhgdMGaZNA

