# The Basics of Python

Python is an interpreted, high-level programming language known for its
flexibility and wide applicability in addition to its accessibility to new
programmers. This is a short description of some key features of Python as they
relate to this project. It is hardly comprehensive, however; if you'd like to
learn more about this language and how to use it, I suggest taking a look at the
[official beginner's guide][1] or this [series of recorded lectures from an MIT
course][2].

## 0 Structure of a Python program
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

## 1 Variables and data types
Nearly all programming languages differentiate between types of data. Different
languages have different standards for how its types are dealt with (Python is
particularly flexible in this regard), but they all require the programmer to
think about, for instance, the difference between the number `12` and the letter
`'a'`. This is partly to help the computer manage all the bits and bytes that
make up a program while it's running, and partly because some operations simply
don't make sense for some types. (What would it mean to compute the sum `12 +
'a'`?)

### 1.1 Numbers
The first type that we'll deal with is the **number**. In Python, "number"
actually comprises four data types that can all mutually interact through
standard arithmetic operations (`+`, `-`, `*` `/`).

#### 1.1.1 `int`
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

#### 1.1.2 `float`
`float`s represent floating-point numbers, or numbers that may have a
non-integer part. `float`s can interact with all other numbers in the way you'd
normally expect, but can be written in two different but functionally equivalent
ways. The first is the standard way, e.g. `250.734`, and the second is with
scientific notation, e.g. `2.50734e2` as shorthand for 2.50734 \* 10^2 or, in
Python syntax, `2.50734 * 10**2`.

#### 1.1.3 `complex`
`complex`s are actually two `float`s bundled together to represent the real and
imaginary parts of a single complex number. To use a `complex`, you can write it
using `j` notation, e.g. `1.5 + 0.4j`, `-2.6 + 1e6j`. When two complex numbers
are added, subtracted, multiplied, or divided together, the imaginary components
obey the usual rules, i.e. `1j * 1j` is equal to `-1`.

#### 1.1.4 `bool`
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

### 1.2 Variables
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

### 1.3 Lists
In the example above, we created a lot of variables with similar names.
`my_num1`, `my_num2`, `my_num3`, and `my_num4` all refer to numbers that are
related to each other, so it would make sense to collect them all together so
that we don't have to keep track of four independent names in the rest of the
program, wouldn't you say? `list`s are the primary way that Python provides to
to this. Specifically, `list`s are ordered, mutable collections of data, all
grouped under one variable name, and are denoted with `[...]` square brackets.
Let's break down what this means.

Python respects a `list`'s order. This means that whenever a list is declared,
for example with
```python
my_list = [1, 5, 3, 8, 5, 0]
```
Python will assign each object contained in the list (in this case the numbers
`1`, `5`, `2`, `8`, `5`, `0`) an index starting with the leftmost at zero. So
under this convention, we'd say that `my_list` contains a `1` at its zeroth
index, a `5` at its first index, and so on:
```python
my_list = [1,  5,  3,  8,  5,  0]
# indices: ^0  ^1  ^2  ^3  ^4  ^5
```
These indices are how Python expects you to refer back to the individual objects
in the list, using the `my_list[...]` notation. To refer to the third index of
`my_list`, you'd write `my_list[3]`. The indices of a `list` can be reassigned
after the `list` is constructed with this notation as well, using the `=` as
above in [Section 1.2](#12-variables): `my_list[3] = 6`.

You can also refer to objects at multiple indices using the `[:]` notation. To
refer to indices `1` through `3` of `my_list`, you'd write `my_list[1:4]`, where
the `:` character (which can only be used in this context) indicates to Python
that you want a range of values, rather than just a single one. Notice that we
have `4` at the upper end of the range, rather than `3`. This is because in
Python, the convention for naming ranges of values is to have the upper end be
non-inclusive; that is, in this case, we name all the indices starting at `1`,
up to *but not including* `4`.

The action of taking only some continuous range of indices from a collection of
objects is called "slicing"; `my_list[1:4]` denotes the "slice" of `my_list`
over indices `1`, `2`, and `3`. In general, you should expect the type of a
slice to be the same as the object that was sliced. In other words, the slice
`my_list[1:4]` should be another list containing the values of `my_list` at
indices `1`, `2`, and `3`.

### 1.4 Strings
`str`s are Python's string type, which represent ordered sequences of characters
that do not hold any numerical value. While the most obvious example of a string
is a collection of alphabetic characters like `"hello world!"`, `str`s may also
contain numbers (as in `"14820"`) but do not behave like numbers in arithmetic
expressions.

`str`s are a lot like `list`s that are denoted with the `"` or `'` symbols
instead of `[` and `]`. On an essential level, they are also an ordered
collection of objects, each of these gets its own index, and the `str` as a
whole follows the `[...]` indexing notation. There are two important
differences, however. The first is that the objects contained in the `str` are
all effectively the same type. Any time you index a `str`, whether by single
index or by slicing, you will get another `str`. The second difference is that
the indices of a `str` are not mutable: it is impossible to alter the value of a
`str` at a specific index through the usual `=` operator.

### 1.5 Dictionaries
`dict`s take the idea of a collection of objects one step further by thinking of
each object as a mapping from an input to an output. Just like you might use a
real dictionary to translate a word in English to one in Spanish, `dict`s take
some "key" and map it to some "value" through the syntax `{key1: value1, key2:
value2, ...}`. The keys and values can be anything you want, but each value in
the `dict` must have at least one corresponding key. Take a look at a couple
examples:
```python
dict_en_sp = {
    "hello": "hola",
    "goodbye": "adiós",
    "my name": "me llamo",
    "where is the library": "donde está la biblioteca"
}

dict_num_word = {
    0: "zero",
    1: "one",
    2: "two",
    3: "three"
}
```



## 2 Logic and control flow

### 2.1 `if` statements

### 2.2 `for` loops

### 2.3 `while` loops

## 3 Functions

## 4 General tips for writing Python

[1]: https://wiki.python.org/moin/BeginnersGuide
[2]: https://www.youtube.com/playlist?list=PLUl4u3cNGP63WbdFxL8giv4yhgdMGaZNA

