Library
=======

All public functions in the ``fastools`` module are directly usable from other
programs. To access the library, simply import ``fastools``:

.. code:: python

    from fastools import gen

    # Make a random FASTA record and write it to `output.fa`.
    handle = open('output.fa', 'w')
    gen(10, handle, 'name', 'description')
