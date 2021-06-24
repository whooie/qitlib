"""
This exercise will be focused mainly on helping you get a feel for working with
basic operations in Python. Here, we'll go over the basic data types of Python,
perform some simple operations and learn how to communicate them in useful ways.

While you're working through this, remember: any result is a useful result, even
if you happen to mess up! There is no 'correct' answer here; the only point of
these exercise is to allow you to play around with this interface.
"""

#== Part 0: `print` ============================================================
print("\n PART 0")
# Let's start with the communication part. The crucial tool for this is the
# `print` function. Just like your everyday mathematical functions, Python
# functions are essentially shorthands for a series of operations that the
# computer performs on some input, and are invoked with
#
# my_function(<my_input>)
#
# `print` is a special function that tells the computer to print out whatever
# its input is. You should use it any time you want to see what your program is
# doing. For example, if you wanted to print out the number 5, you'd wrap it in
# parentheses and call `print` on it, like so:
#
# print(5)
#
# So let's try a few  examples. Below, we have a few things we'd # like to print
# out. But if we were to just run this program as-is, we wouldn't # see anything
# because we're not telling the computer exactly what we want. To # properly let
# it know, we need to use the `print` function. Try calling the `print` function
# on each of these things.

"hello world!"

3.1415926535897932384626

["Alice", "Bob", "Charlie", "Debora", "Emily"]

# Great! Now try running this by using the `cd` command in the console to get to
# this folder and entering `python 00-basic_python.py`. You should be able to
# see whatever's inside the parentheses pop up in the panel on the right.
#
# Now try printing some other things out. They can be whatever you want as long
# as they are numbers or words surrounded by quotes, like "hello world!" above.

# --- your code here: ---





# -----------------------



#== PART 1: Numbers and math ===================================================
print("\nPART 1")
# Python has a few different data types we'll be working with. Most of them are
# numbers, so we'll start there with some basic arithmetic. First, let's verify
# that Python can actually do math. How should we go about this? Well, there are
# four basic arithmetic operations (addition, subtraction, multiplication, and
# division), which are indicated by the `+`, `-`, `*`, and `/` characters, so we
# should begin by testing out Python's ability to perform them. Try printing out
# the following operations, and make sure they agree with what you expect:

1 + 1

6 - 8.2

8 * 5

5 / 2.5

# We can also test Python's ability to follow the correct order of operations
# with more complicated expressions. Try printing these out and making sure they
# also agree with what you get when you work them out by hand:

5 * 2 - 10 / 18

(7 + 8) / 10 * (1 + 10)

5 * (5 + 6 / 3 - 10)

# One other thing we can use to test out Python math is division by zero. We can
# divide almost any number we like, but what happens when we try dividing by
# zero? Try printing out the result of `1 / 0`:

# --- your code here: ---



# -----------------------

# If you look at the output of your program, you should be seeing something like
# `ZeroDivisionError: division by zero`. When Python encounters operations like
# this that generate errors, it immediately halts the execution of the program
# because it doesn't know what to do. So the system works! Except now we have to
# deal with the offending operation, otherwise nothing past that point would
# ever get executed. One solution is to simply delete that line, but another
# thing you could do is to use "comments". Python has a built-in functionality
# where anything on a line following a `#` character is completely ignored
# (which is why you can read all this text without it affecting your program).
# Programmers often use this to leave notes in their code so that other
# programmers can quickly understand what the code is doing, but it can also be
# a tool to "disable" certain lines in the program so that errors don't pop up.
# Try adding a `#` right at the beginning of your code printing out the result
# of `1 / 0` and re-running the program. If you did it correctly, everything
# should run without the `ZeroDivisionError` popping up.



#== PART 2: Variables ==========================================================
print("\nPART 2")
# So you've learned how to do math and print out the results. Great! But what if
# you wanted to store your results for later use? Let's say you wanted to print
# out successive powers of 5. One way to do this is through the simple, direct
# way that we've been looking at so far:
#
# print(5)
# print(5 * 5)
# print(5 * 5 * 5)
# print(5 * 5 * 5 * 5)
# ...
#
# and so on. But that's a lot of wasted computation. Since we know that any
# power of 5 is just the previous power of 5 multiplied by another 5, we're
# making the computer do unnecessary work when we tell it to multiply all those
# 5's together for each `print` statement. A much better way to do this is with
# variables.
#
# Simply put, a variable (in Python) is a name you can give any value so that
# you can refer back to it at a later time. Assigning values to variables is
# done through, unsurprisingly, the "variable assignment" operator, `=`. The
# general form of this operation looks like this:
#
# my_variable = <value>
#
# where <value> can be anything you want, as long as it's a valid object in
# Python. Try setting some of your own variables below:

# --- your code here ---





# ----------------------

# Once you have a value assigned to a variable, you can use the variable just
# like you would use its value. For example:
#
# my_number = 5
# print(2 * my_number + 1)
#
# Here, `my_number` stands in place for whatever value we've assigned to it,
# which in this case is 5. This means we can even use variables to assign values
# to other variables:
#
# my_number = 5
# my_other_number = 2 * my_number + 1
# print(my_other_number)
#
# This is equivalent to the other example above! You can also assign a variable
# to the result of an expression involving itself:
#
# my_number = 1
# my_number = my_number + 2
#
# On the second line, we're taking whatever the value of `my_number` is at that
# time (in this case, `my_number` is 1), adding 2 to it, and then using it to
# overwrite the previous value assigned to `my_number`. (What is this value?)
#
# As another example, try thinking about what will get printed out by the
# following:

#my_number = 2
#print(my_number)
#my_number = my_number + 2
#print(my_number)
#my_number = my_number + 2
#print(my_number)

# Once you've thought about it, try deleting the `#` characters and seeing if
# you were right! Finally, try using variables to find a better way to print out
# successive powers of 5 than the one at the beginning of this section. Hint: it
# will look a lot like what you've just uncommented above.

# --- your code here ---





# ----------------------



#== PART 3: Containers =========================================================
print("\nPART 3")
# In Python, there's a broad category of things called "containers" --
# essentially, objects that contain other objects by acting like collections of
# variables. A few of the most commonly used containers are `list`s, `str`s, and
# `dict`s. Let's start with the first of these.
#
# A `list` is, in technical terms, an ordered, mutable collection of data. Let's
# break down what this means. To be a "collection of data" is fairly
# straightforward: `list`s hold a collection of other pieces of data. These
# other data objects can be anything we want, and the `list` functions basically
# as a way to group them all under one variable name, with a way to access each
# object individually or in sub-groups.
#
# Being "ordered" means that Python respects the order in which the objects were
# placed in the `list`. Whenever a list is declared, for example with

my_list = [1, 5, 3, 8, 5, 0]

# (notice the square brackets and the commas separating each number), Python
# will assign each object contained in the `list` an index starting with the
# leftmost one at zero. So under this convention, we'd say that `my_list`
# contains a 1 at its zeroth index, a 5 at its first index, and so on:

my_list = [1,  5,  3,  8,  5,  0]
# indices: ^0  ^1  ^2  ^3  ^4  ^5

# These indices are how Python expects you to refer back to the individual
# objects in the `list`, using the `my_list[...]` notation. To refer to the
# third index of `my_list`, you'd write `my_list[3]`. You can also refer to
# multiple indices at once using the `[:]` notation. To refer to indices `1`
# through `3` of `my_list`, you'd write `my_list[1:4]`, where the `:` character
# indicates to Python that you want a range of values, rather than just a single
# one. Notice that we have `4` at the upper end of the range, rather than `3`.
# This is because in Python, the convection for naming ranges of values is to
# have the upper end be non-inclusive: we name all the indices starting at `1`
# and going up to *but not including* `4`.
#
# Try wrapping the following in `print` calls. For each one, try to think about
# what you should be seeing before you hit "run".

my_list[3]
my_list[1]
my_list[5]
my_list[1:4]
my_list[0:6]
my_list[4:5]

# Additionally, being "mutable" means that both the objects held by the `list`
# and the `list` itself can be altered after it's created. We won't be doing
# much with this, but I encourage you to look into how it's accomplished if
# you're interested!
#
# Next up are `str`s. "str" is short for "string", which represent ordered
# sequences of characters that do not hold any numerical value. Things like

my_str1 = "hello world!"
my_str2 = "!$@#%"
my_str3 = '3.1415926535897932384626'
my_str4 = '你好世界'

# are all strings, indicated by being wrapped in either double quotes (`"`) or
# single quotes(`'`). As you can see, `str`s may also contain numbers (as in the
# third string listed here), but will not behave as numbers -- for example, you
# can't do this:
#
# "4" + 5
#
# `str`s are a lot like `list`s. On an essential level they are also an ordered
# collection of objects, so you can use the same `[...]` indexing notation as
# for `list`s on them. But there are two important differences. The first is
# that all the objects contained in a `str` are all the same, and consist of
# only a single character. Any time you index a `str`, whether by a single
# index or with a range of indices, you will get another `str` out. The second
# difference is that `str`s are not mutable.
#
# Try performing some indexing with `my_str1`, `my_str2`, `my_str3`, and
# `my_str4` defined above.

# --- your code here ---





# ----------------------

# Finally, we have `dict`, which is short for "dictionary". Just like you might
# use a real dictionary to translate a word in English to one in Spanish,
# `dict`s some "key" (the English word) and relate it to some "value" (the
# Spanish word). The way to write a `dict` is more complicated than for `list`s
# and `str`s, but it looks like this:
#
# {key1: value1, key2: value2, key3: value3, ...}
#
# Notice that each key is paired with its value through a `:` character, and
# separated from all the other key-value pairs with a `,`, while the whole thing
# is wrapped in `{` / `}` curly braces. Let's take a look at a couple examples.

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

# The first dictionary above pairs `str`s together with `str`s, and the second
# pairs numbers with `str`s. To make use of these pairings, we use the `[...]`
# indexing syntax from before, except we can feed in keys (which can be
# anything) instead of only indices. Try printing out the following:

dict_en_sp["hello"]
dict_en_sp["goodbye"]
dict_num_word[0]
dict_num_word[3]

# The important difference with this `[...]` notation is that `dict`s cannot
# accept ranges of keys like `list`s can. Why? Because the other important
# difference between `list`s and `dict`s is that `dict`s are *not* ordered. A
# good way to think about it is that `list`s are almost like `dict`s where all
# the keys have to be numbers, beginning with zero, which naturally places
# everything in an order. But since the keys can be anything for `dict`s,
# there's no way to apply an ordering, so `dict`s in fact *cannot* be ordered.
#
# On the other hand, though, `dict`s are mutable like `list`s. To add more
# key-value pairings, all you have to do is assign them like variables:
#
# my_dict[<my_new_key>] = <my_new_value>
#
# So for `dict_en_sp` and `dict_num_word`, we could add a couple more pairings
# like so:

dict_en_sp["the spider"] = "la araña"
dict_en_sp["the python"] = "la pitón"
dict_num_word[4] = "four"
dict_num_word[100] = "one hundred"

# Try adding a few of your own, and then call `print` on both dictionaries to
# see what happens!

# --- your code here ---





#print(dict_en_sp)
#print(dict_num_word)
# ----------------------



#== PART 4: Functions and methods ==============================================
print("\nPART 4")
# Functions are an integral part of any kind of useful programming.
