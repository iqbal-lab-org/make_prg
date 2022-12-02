This file describes this project's coding conventions.

1. Do not directly commit to dev or main/master. Do a PR with at least one reviewer;

2. Every code should be covered by at least one test, except for justifiable exceptions (e.g. argument parsing code);

3. Asserts should be used to check/avoid developer errors, while exceptions should be used for errors that might
appear due to incorrect usage of the tool by the user (e.g. a badly formed fasta file). In assert statements, please
add a detailed description explaining why we think this particular event should not happen, unless the assertion message
or condition tested is already clear enough.

4. TODO: add lots more.
