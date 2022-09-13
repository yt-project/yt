A Brief Introduction to Python
------------------------------

All scripts that use yt are really Python scripts that use yt as a library.
The great thing about Python is that the standard set of libraries that come
with it are very extensive -- Python comes with everything you need to write
and run mail servers and web servers, create Logo-style turtle graphics, do
arbitrary precision math, interact with the operating system, and many other
things.  In addition to that, efforts by the scientific community to improve
Python for computational science have created libraries for fast array
computation, GPGPU operations, distributed computing, and visualization.

So when you use yt through the scripting interface, you get for free the
ability to interlink it with any of the other libraries available for Python.
In the past, this has been used to create new types of visualization using
OpenGL, data management using Google Docs, and even a simulation that sends an
SMS when it has new data to report on.

But, this also means learning a little bit of Python!  This next section
presents a short tutorial of how to start up and think about Python, and then
moves on to how to use Python with yt.

Starting Python
+++++++++++++++

Python has two different modes of execution: interactive execution, and
scripted execution.  We'll start with interactive execution and then move on to
how to write and use scripts.

Before we get started, we should briefly touch upon the commands ``help`` and
``dir``.  These two commands provide a level of introspection:
``help(something)`` will return the internal documentation on ``something``,
including how it can be used and all of the possible "methods" that can be
called on it.  ``dir()`` will return the available commands and objects that
can be directly called, and ``dir(something)`` will return information about
all the commands that ``something`` provides.  This probably sounds a bit
opaque, but it will become clearer with time -- it's also probably helpful to
call ``help`` on any or all of the objects we create during this orientation.

To start up Python, at your prompt simply type:

.. code-block:: bash

  $ python

This will open up Python and give to you a simple prompt of three greater-than
signs.  Let's inaugurate the occasion appropriately -- type this::

   >>> print("Hello, world.")

As you can see, this printed out the string "Hello, world." just as we
expected.  Now let's try a more advanced string, one with a number in it.  For
this we'll use an "f-string", which is the preferred way to format strings in modern Python.
We'll print pi, but only with three digits of accuracy.::

   >>> print(f"Pi is precisely {3.1415926:0.2f}")

This took the number we fed it (3.1415926) and printed it out as a floating
point number with two decimal places. Now let's try something a bit different
-- let's print out both the name of the number and its value.::

   >>> print(f"{'pi'} is precisely {3.1415926:0.2f}")

And there you have it -- the very basics of starting up Python, and some very
simple mechanisms for printing values out.  Now let's explore a few types of
data that Python can store and manipulate.

Data Types
++++++++++

Python provides a number of datatypes, but the main ones that we'll concern
ourselves with at first are lists, tuples, strings, numbers, and dictionaries.
Most of these can be instantiated in a couple different ways, and we'll look at
a few of them.  Some of these objects can be modified in place, which is called
being mutable, and some are immutable and cannot be modified in place.  We'll
talk below about what that means.

Perhaps most importantly, though, is an idea about how Python works in terms of
names and bindings -- this is called the "object model."  When you create an
object, that is independent from binding it to a name -- think of this like
pointers in C.  This also operates a bit differently for mutable and immutable
types.  We'll talk a bit more about this later, but it's handy to initially
think of things in terms of references.  (This is also, not coincidentally, how
Python thinks of things internally as well!)  When you create an object,
initially it has no references to it -- there's nothing that points to it.
When you bind it to a name, then you are making a reference, and so its
"reference count" is now 1.  When something else makes a reference to it, the
reference count goes up to 2.  When the reference count returns to 0, the
object is deleted -- but not before.  This concept of reference counting comes
up from time to time, but it's not going to be a focus of this orientation.

The two easiest datatypes are simply strings and numbers.  We can make a string
very easily::

   >>> my_string = "Hello there"
   >>> print(my_string)

We can also take a look at each individual part of a string.  We'll use the
'slicing' notation for this.  As a brief note, slicing is 0-indexed, so that
element 0 corresponds to the first element.  If we wanted to see the third
element of our string::

   >>> print(my_string[2])

We can also take the third through the 5 elements::

   >>> print(my_string[2:5])

But note that if you try to change an element directly, Python objects and it
won't let you -- that's because strings are immutable.  (But, note that because
of how the += operator works, we can do "my_string += '1'" without issue.)

To create a number, we do something similar::

   >>> a = 10
   >>> print(a)

This works for floating points as well.  Now we can do math on these numbers::

   >>> print(a**2)
   >>> print(a + 5)
   >>> print(a + 5.1)
   >>> print(a / 2.0)


Now that we have a couple primitive datatypes, we can move on to sequences --
lists and tuples.  These two objects are very similar, in that they are
collections of arbitrary data types.  We'll only look at collections of strings
and numbers for now, but these can be filled with arbitrary datatypes
(including objects that yt provides, like spheres, datasets, grids, and
so on.)  The easiest way to create a list is to simply construct one::

   >>> my_list = []

At this point, you can find out how long it is, you can append elements, and
you can access them at will::

   >>> my_list.append(1)
   >>> my_list.append(my_string)
   >>> print(my_list[0])
   >>> print(my_list[-1])
   >>> print(len(my_list))

You can also create a list already containing an initial set of elements::

   >>> my_list = [1, 2, 3, "four"]
   >>> my_list[2] = "three!!"

Lists are very powerful objects, which we'll talk about a bit below when
discussing how iteration works in Python.

A tuple is like a list, in that it's a sequence of objects, and it can be
sliced and examined piece by piece.  But unlike a list, it's immutable:
whatever a tuple contains at instantiation is what it contains for the rest of
its existence.  Creating a tuple is just like creating a list, except that you
use parentheses instead of brackets::

   >>> my_tuple = (1, "a", 62.6)

Tuples show up very commonly when handling arguments to Python functions and
when dealing with multiple return values from a function.  They can also be
unpacked::

   >>> v1, v2, v3 = my_tuple

will assign 1, "a", and 62.6 to v1, v2, and v3, respectively.

Mutables vs Immutables and Is Versus Equals
+++++++++++++++++++++++++++++++++++++++++++

This section is not a "must read" -- it's more of an exploration of how
Python's objects work.  At some point this is something you may want to be
familiar with, but it's not strictly necessary on your first pass.

Python provides the operator ``is`` as well as the comparison operator ``==``.
The operator ``is`` determines whether two objects are in fact the same object,
whereas the operator ``==`` determines if they are equal, according to some
arbitrarily defined equality operation.  Think of this like comparing the
serial numbers on two pictures of a dollar bill (the ``is`` operator) versus
comparing the values of two pieces of currency (the ``==`` operator).

This digs in to the idea of how the Python object model works, so let's test
some things out.  For instance, let's take a look at comparing two floating
point numbers::

   >>> a = 10.1
   >>> b = 10.1
   >>> print(a == b)
   >>> print(a is b)

The first one returned True, but the second one returned False.  Even though
both numbers are equal, they point to different points in memory.  Now let's
try assigning things a bit differently::

   >>> b = a
   >>> print(a is b)

This time it's true -- they point to the same part of memory.  Try incrementing
one and seeing what happens.  Now let's try this with a string::

   >>> a = "Hi there"
   >>> b = a
   >>> print(a is b)

Okay, so our intuition here works the same way, and it returns True.  But what
happens if we modify the string?::

   >>> a += "!"
   >>> print(a)
   >>> print(b)
   >>> print(a is b)

As you can see, now not only does a contain the value "Hi there!", but it also
is a different value than what b contains, and it also points to a different
region in memory.  That's because strings are immutable -- the act of adding on
"!" actually creates an entirely new string and assigns that entirely new
string to the variable a, leaving the string pointed to by b untouched.

With lists, which are mutable, we have a bit more liberty with how we modify
the items and how that modifies the object and its pointers.  A list is really
just a pointer to a collection; the list object itself does not have any
special knowledge of what constitutes that list.  So when we initialize a and
b::

   >>> a = [1, 5, 1094.154]
   >>> b = a

We end up with two pointers to the same set of objects.  (We can also have a
list inside a list, which adds another fun layer.)  Now when we modify a, it
shows up in b::

   >>> a.append("hat wobble")
   >>> print(b[-1])

This also works with the concatenation operator::

   >>> a += ["beta sequences"]
   >>> print(a[-1], b[-1])

But we can force a break in this by slicing the list when we initialize::

   >>> a = [1, 2, 3, 4]
   >>> b = a[:]
   >>> a.append(5)
   >>> print(b[-1], a[-1])

Here they are different, because we have sliced the list when initializing b.

The coolest datatype available in Python, however, is the dictionary.  This is
a mapping object of key:value pairs, where one value is used to look up another
value.  We can instantiate a dictionary in a variety of ways, but for now we'll
only look at one of the simplest mechanisms for doing so::

   >>> my_dict = {}
   >>> my_dict["A"] = 1.0
   >>> my_dict["B"] = 154.014
   >>> my_dict[14001] = "This number is great"
   >>> print(my_dict["A"])

As you can see, one value can be used to look up another.  Almost all datatypes
(with a few notable exceptions, but for the most part these are quite uncommon)
can be used as a key, and you can use any object as a value.

We won't spend too much time discussing dictionaries explicitly, but I will
leave you with a word on their efficiency: the Python lookup algorithm is known
for its hand-tuned optimization and speed, and it's very common to use
dictionaries to look up hundreds of thousands or even millions of elements and
to expect it to be responsive.

Looping
+++++++

Looping in Python is both different and more powerful than in lower-level
languages.  Rather than looping based exclusively on conditionals (which is
possible in Python) the fundamental mode of looping in Python is iterating
over objects.  In C, one might construct a loop where some counter variable is
initialized, and at each iteration of the loop it is incremented and compared
against a reference value; when the counter variable reaches the reference
variable, the loop is terminated.

In Python, on the other hand, to accomplish iteration through a set of
sequential integers, one actually constructs a sequence of those integers, and
iterates over that sequence.  For more discussion of this, and some very, very
powerful ways of accomplishing this iteration process, look through the Python
documentation for the words 'iterable' and 'generator.'

To see this in action, let's first take a look at the built-in function
``range``. ::

   >>> print(range(10))

As you can see, what the function ``range`` returns is a list of integers,
starting at zero, that is as long as the argument to the ``range`` function.
In practice, this means that calling ``range(N)`` returns ``0, 1, 2, ... N-1``
in a list.  So now we can execute a for loop, but first, an important
interlude:

Control blocks in Python are delimited by white space.

This means that, unlike in C with its brackets, you indicate an isolated
control block for conditionals, function declarations, loops and other things
with an indentation.  When that control block ends, you dedent the text.  In
yt, we use four spaces -- I recommend you do the same -- which can be inserted
by a text editor in place of tab characters.

Let's try this out with a for loop.  First type ``for i in range(10):`` and
press enter.  This will change the prompt to be three periods, instead of three
greater-than signs, and you will be expected to hit the tab key to indent.
Then type "print(i)", press enter, and then instead of indenting again, press
enter again.  The entire entry should look like this::

   >>> for i in range(10):
   ...     print(i)
   ...

As you can see, it prints out each integer in turn.  So far this feels a lot
like C.  (It won't, if you start using iterables in place of sequences -- for
instance, ``xrange`` operates just like range, except instead of returning an
already-created list, it returns the promise of a sequence, whose elements
aren't created until they are requested.)  Let's try it with our earlier list::

   >>> my_sequence = ["a", "b", 4, 110.4]
   >>> for i in my_sequence:
   ...     print(i)
   ...

This time it prints out every item in the sequence.

A common idiom that gets used a lot is to figure out which index the loop is
at.  The first time this is written, it usually goes something like this::

   >>> index = 0
   >>> my_sequence = ["a", "b", 4, 110.4]
   >>> for i in my_sequence:
   ...     print("%s = %s" % (index, i))
   ...     index += 1
   ...

This does what you would expect: it prints out the index we're at, then the
value of that index in the list.  But there's an easier way to do this, less
prone to error -- and a bit cleaner!  You can use the ``enumerate`` function to
accomplish this::

   >>> my_sequence = ["a", "b", 4, 110.4]
   >>> for index, val in enumerate(my_sequence):
   ...     print("%s = %s" % (index, val))
   ...

This does the exact same thing, but we didn't have to keep track of the counter
variable ourselves.  You can use the function ``reversed`` to reverse a
sequence in a similar fashion.  Try this out::

   >>> my_sequence = range(10)
   >>> for val in reversed(my_sequence):
   ...     print(val)
   ...

We can even combine the two!::

   >>> my_sequence = range(10)
   >>> for index, val in enumerate(reversed(my_sequence)):
   ...     print("%s = %s" % (index, val))
   ...

The most fun of all the built-in functions that operate on iterables, however,
is the ``zip`` function.  This function will combine two sequences (but only up
to the shorter of the two -- so if one has 16 elements and the other 1000, the
zipped sequence will only have 16) and produce iterators over both.

As an example, let's say you have two sequences of values, and you want to
produce a single combined sequence from them.::

   >>> seq1 = ["Hello", "What's up", "I'm fine"]
   >>> seq2 = ["!", "?", "."]
   >>> seq3 = []
   >>> for v1, v2 in zip(seq1, seq2):
   ...     seq3.append(v1 + v2)
   ...
   >>> print(seq3)

As you can see, this is much easier than constructing index values by hand and
then drawing from the two sequences using those index values.  I should note
that while this is great in some instances, for numeric operations, NumPy
arrays (discussed below) will invariably be faster.

Conditionals
++++++++++++

Conditionals, like loops, are delimited by indentation.  They follow a
relatively simple structure, with an "if" statement, followed by the
conditional itself, and then a block of indented text to be executed in the
event of the success of that conditional.  For subsequent conditionals, the
word "elif" is used, and for the default, the word "else" is used.

As a brief aside, the case/switch statement in Python is typically executed
using an if/elif/else block; this can be done using more complicated
dictionary-type statements with functions, but that typically only adds
unnecessary complexity.

For a simple example of how to do an if/else statement, we'll return to the
idea of iterating over a loop of numbers.  We'll use the ``%`` operator, which
is a binary modulus operation: it divides the first number by the second and
then returns the remainder.  Our first pass will examine the remainders from
dividing by 2, and print out all the even numbers.  (There are of course easier
ways of determining which numbers are multiples of 2 -- particularly using
NumPy, as we'll do below.)::

   >>> for val in range(100):
   ...     if val % 2 == 0:
   ...         print("%s is a multiple of 2" % (val))
   ...

Now we'll add on an ``else`` statement, so that we print out all the odd
numbers as well, with the caveat that they are not multiples of 2.::

   >>> for val in range(100):
   ...     if val % 2 == 0:
   ...         print("%s is a multiple of 2" % (val))
   ...     else:
   ...         print("%s is not a multiple of 2" % (val))
   ...

Let's extend this to check the remainders of division with both 2 and 3, and
determine which numbers are multiples of 2, 3, or neither.  We'll do this for
all numbers between 0 and 99.::

   >>> for val in range(100):
   ...     if val % 2 == 0:
   ...         print("%s is a multiple of 2" % (val))
   ...     elif val % 3 == 0:
   ...         print("%s is a multiple of 3" % (val))
   ...     else:
   ...         print("%s is not a multiple of 2 or 3" % (val))
   ...

This should print out which numbers are multiples of 2 or 3 -- but note that
we're not catching all the multiples of 6, which are multiples of both 2 and 3.
To do that, we have a couple options, but we can start with just changing the
first if statement to encompass both, using the ``and`` operator::

   >>> for val in range(100):
   ...     if val % 2 == 3 and val % 3 == 0:
   ...         print("%s is a multiple of 6" % (val))
   ...     elif val % 2 == 0:
   ...         print("%s is a multiple of 2" % (val))
   ...     elif val % 3 == 0:
   ...         print("%s is a multiple of 3" % (val))
   ...     else:
   ...         print("%s is not a multiple of 2 or 3" % (val))
   ...

In addition to the ``and`` statement, the ``or`` and ``not`` statements work in
the expected manner.  There are also several built-in operators, including
``any`` and ``all`` that operate on sequences of conditionals, but those are
perhaps better saved for later.

Array Operations
++++++++++++++++

In general, iteration over sequences carries with it some substantial overhead:
each value is selected, bound to a local name, and then its type is determined
when it is acted upon.  This is, regrettably, the price of the generality that
Python brings with it.  While this overhead is minimal for operations acting on
a handful of values, if you have a million floating point elements in a
sequence and you want to simply add 1.2 to all of them, or multiply them by
2.5, or exponentiate them, this carries with it a substantial performance hit.

To accommodate this, the NumPy library has been created to provide very fast
operations on arrays of numerical elements.  When you create a NumPy array, you
are creating a shaped array of (potentially) sequential locations in memory
which can be operated on at the C-level, rather than at the interpreted Python
level.  For this reason, which NumPy arrays can act like Python sequences can,
and can thus be iterated over, modified in place, and sliced, they can also be
addressed as a monolithic block.  All of the fluid and particle quantities used
in yt will be expressed as NumPy arrays, allowing for both efficient
computation and a minimal memory footprint.

For instance, the following operation will not work in standard Python::

   >>> vals = range(10)
   >>> vals *= 2.0

(Note that multiplying vals by the integer 2 will not do what you think: rather
than multiplying each value by 2.0, it will simply double the length of the
sequence!)

To get started with array operations, let's first import the NumPy library.
This is the first time we've seen an import in this orientation, so we'll
dwell for a moment on what this means.  When a library is imported, it is read
from disk, the functions are loaded into memory, and they are made available
to the user.  So when we execute::

   >>> import numpy

The ``numpy`` module is loaded, and then can be accessed::

   >>> numpy.arange(10)

This calls the ``arange`` function that belongs to the ``numpy`` module's
"namespace."  We'll use the term namespace to refer to the variables,
functions, and submodules that belong to a given conceptual region.  We can
also extend our current namespace with the contents of the ``numpy`` module, so
that we don't have to prefix all of our calling of ``numpy`` functions with
``numpy.`` but we will not do so here, so as to preserve the distinction
between the built-in Python functions and the NumPy-provided functions.

To get started, let's perform the NumPy version of getting a sequence of
numbers from 0 to 99::

   >>> my_array = numpy.arange(100)
   >>> print(my_array)
   >>> print(my_array * 2.0)
   >>> print( my_array * 2)

As you can see, each of these operations does exactly what we think it ought
to.  And, in fact, so does this one::

   >>> my_array *= 2.0

So far we've only examined what happens if we have operate on a single array of
a given shape -- specifically, if we have an array that is N elements long, but
only one dimensional.  NumPy arrays are, for the most part, defined by their
data, their shape, and their data type.  We can examine both the shape (which
includes dimensionality) and the size (strictly the total number of elements)
in an array by looking at a couple properties of the array::

   >>> print(my_array.size)
   >>> print(my_array.shape)

Note that size must be the product of the components of the shape.  In this
case, both are 100.  We can obtain a new array of a different shape by calling
the ``reshape`` method on an array::

   >>> print(my_array.reshape((10, 10)))

In this case, we have not modified ``my_array`` but instead created a new array
containing the same elements, but with a different dimensionality and shape.
You can modify an array's shape in place, as well, but that should be done with
care and the explanation of how that works and its caveats can come a bit
later.

There are a few other important characteristics of arrays, and ways to create
them.  We can see what kind of datatype an array is by examining its ``dtype``
attribute::

   >>> print(my_array.dtype)

This can be changed by calling ``astype`` with another datatype.  Datatypes
include, but are not limited to, ``int32``, ``int64``, ``float32``,
``float64``.::

   >>> float_array = my_array.astype("float64")

Arrays can also be operated on together, in lieu of something like an iteration
using the ``zip`` function.  To show this, we'll use the
``numpy.random.random`` function to generate a random set of values of length
100, and then we'll multiply our original array against those random values.::

   >>> rand_array = numpy.random.random(100)
   >>> print(rand_array * my_array)

There are a number of functions you can call on arrays, as well.  For
instance::

   >>> print(rand_array.sum())
   >>> print(rand_array.mean())
   >>> print(rand_array.min())
   >>> print(rand_array.max())

Indexing in NumPy is very fun, and also provides some advanced functionality
for selecting values.  You can slice and dice arrays::

   >>> print(my_array[50:60])
   >>> print(my_array[::2])
   >>> print(my_array[:-10])

But Numpy also provides the ability to construct boolean arrays, which are the
result of conditionals.  For example, let's say that you wanted to generate a
random set of values, and select only those less than 0.2::

   >>> rand_array = numpy.random.random(100)
   >>> print(rand_array < 0.2)

What is returned is a long list of booleans.  Boolean arrays can be used as
indices -- what this means is that you can construct an index array and then
use that toe select only those values where that index array is true.  In this
example we also use the ``numpy.all`` and ``numpy.any`` functions, which do
exactly what you might think -- they evaluate a statement and see if all
elements satisfy it, and if any individual element satisfies it,
respectively.::

   >>> ind_array = rand_array < 0.2
   >>> print(rand_array[ind_array])
   >>> print(numpy.all(rand_array[ind_array] < 0.2))

You can even skip the creation of the variable ``ind_array`` completely, and
instead just coalesce the statements into a single statement::

   >>> print(numpy.all(rand_array[rand_array < 0.2] < 0.2))
   >>> print(numpy.any(rand_array[rand_array < 0.2] > 0.2))

You might look at these and wonder why this is useful -- we've already selected
those elements that are less than 0.2, so why do we want to re-evaluate it?
But the interesting component to this is that a conditional applied to one
array can be used to index another array.  For instance::

   >>> print(my_array[rand_array < 0.2])

Here we've identified those elements in our random number array that are less
than 0.2, and printed the corresponding elements from our original sequential
array of integers.  This is actually a great way of selecting a random sample
of a dataset -- in this case we get back approximately 20% of the dataset
``my_array``, selected at random.

To create arrays from nothing, several options are available.  The command
``numpy.array`` will create an array from any arbitrary sequence::

   >>> my_sequence = [1.0, 510.42, 1789532.01482]
   >>> my_array = numpy.array(my_sequence)

Additionally, arrays full of ones and zeros can be created::

   >>> my_integer_ones = numpy.ones(100)
   >>> my_float_ones = numpy.ones(100, dtype="float64")
   >>> my_integer_zeros = numpy.zeros(100)
   >>> my_float_zeros = numpy.zeros(100, dtype="float64")

The function ``numpy.concatenate`` is also useful, but outside the scope of
this orientation.

The NumPy documentation has a number of more advanced mechanisms for combining
arrays; the documentation for "broadcasting" in particular is very useful, and
covers mechanisms for combining arrays of different shapes and sizes, which can
be tricky but also extremely powerful.  We won't discuss the idea of
broadcasting here, simply because I don't know that I could do it justice!  The
NumPy Docs have a great `section on broadcasting
<https://numpy.org/doc/stable/user/basics.broadcasting.html>`_.

Scripted Usage
++++++++++++++

We've now explored Python interactively.  However, for long-running analysis
tasks or analysis tasks meant to be run on a compute cluster non-interactively,
we will want to utilize its scripting interface.  Let's start by quitting out
of the interpreter.  If you have not already done so, you can quit by pressing
"Ctrl-D", which will free all memory used by Python and return you to your
shell's command prompt.

At this point, open up a text editor and edit a file called
``my_first_script.py``.  Python scripts typically end in the extension ``.py``.
We'll start our scripting tests by doing some timing of array operations versus
sequence operations.  Into this file, type this text::

   import numpy
   import time

   my_array = numpy.arange(1000000, dtype="float64")

   t1 = time.time()
   my_array_squared = my_array**2.0
   t2 = time.time()

   print("It took me %0.3e seconds to square the array using NumPy" % (t2-t1))

   t1 = time.time()
   my_sequence_squared = []
   for i in range(1000000):
       my_sequence_squared.append(i**2.0)
   t2 = time.time()

   print("It took me %0.3e seconds to square the sequence without NumPy" % (t2-t1))

Now save this file, and return to the command prompt.  We can execute it by
supplying it to Python:

.. code-block:: bash

   $ python my_first_script.py

It should run, display two pieces of information, and terminate, leaving you
back at the command prompt.  On my laptop, the array operation is approximately
42 times faster than the sequence operation!  Of course, depending on the
operation conducted, this number can go up quite substantially.

If you want to run a Python script and then be given a Python interpreter
prompt, you can call the ``python`` command with the option ``-i``:

.. code-block:: bash

   $ python -i my_first_script.py

Python will execute the script and when it has reached the end it will give you
a command prompt.  At this point, all of the variables you have set up and
created will be available to you -- so you can, for instance, print out the
contents of ``my_array_squared``::

   >>> print(my_array_squared)

The scripting interface for Python is quite powerful, and by combining it with
interactive execution, you can, for instance, set up variables and functions
for interactive exploration of data.

Functions and Objects
+++++++++++++++++++++

Functions and Objects in Python are the easiest way to perform very complex,
powerful actions in Python.  For the most part we will not discuss them; in
fact, the standard Python tutorial that comes with the Python documentation is
a very good explanation of how to create and use objects and functions, and
attempting to replicate it here would simply be futile.

yt provides both many objects and functions for your usage, and it is through
the usage and combination of functions and objects that you will be able to
create plots, manipulate data, and visualize your data.

And with that, we conclude our brief introduction to Python.  I recommend
checking out the standard Python tutorial or browsing some of the NumPy
documentation.  If you're looking for a book to buy, the only book I've
personally ever been completely satisfied with has been David Beazley's book on
Python Essentials and the Python standard library, but I've also heard good
things about many of the others, including those by Alex Martelli and Wesley
Chun.

We'll now move on to talking more about how to use yt, both from a scripting
perspective and interactively.

Python and Related References
+++++++++++++++++++++++++++++
    * `Python quickstart <https://docs.python.org/3/tutorial/>`_
    * `Learn Python the Hard Way <https://learnpythonthehardway.org/python3/>`_
    * `Byte of Python <https://python.swaroopch.com/>`_
    * `Dive Into Python <https://diveintopython3.problemsolving.io/>`_
    * `Numpy docs <https://numpy.org/doc/stable/>`_
    * `Matplotlib docs <https://matplotlib.org>`_
