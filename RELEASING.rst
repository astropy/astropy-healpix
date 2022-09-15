Release instructions
====================

Once the package is ready to release, use ``git tag`` to tag the
release::

    git tag -m <version> <version>

e.g::

    git tag -m v0.1 v0.1

You can also include the ``-s`` flag to sign the tag if you have
PGP keys set up. Then, push the tag to GitHub, e.g.::

    git push upstream v0.1

and the build should happen automatically on GitHub Actions. You can
follow the progress of the build here:

https://github.com/astropy/astropy-healpix/actions/workflows/publish.yml

If there are any failures, you can always delete the tag, fix the
issues, tag the release again, and push the tag to GitHub.

See the `OpenAstronomy GitHub Actions Workflows
<https://github.com/OpenAstronomy/github-actions-workflows>`_
for more details about the GitHub Actions set-up.
