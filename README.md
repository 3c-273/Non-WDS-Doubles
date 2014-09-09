Non-WDS-Doubles
===============

Find Doubles not in the WDS

Here are C and Perl programs that will:

 -> Make a list of all WDS entry's coordinates and list them in order of increasing RA.

 -> Parse all entries in the UCAC4 and split them into files that each contain about a square degree of sky

 -> Create a list of candidates stars from the UCAC4, all brighter than a specified magnitude.

 -> Search each candidate to see if stars near it could be possible secondaries, based on separation, brightness, brightness differential between the primary and secondary, amount of proper motion and similar proper motion.

 -> Check to see if a candidate pairs is already in the WDS.

 -> Make an HTML formatted  list of candidate pairs that are not in the WDS and have passed all of the above criteria.
