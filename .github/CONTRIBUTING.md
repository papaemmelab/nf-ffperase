# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given. This project could always use more documentation, whether as part of the README, in docstrings, or even on the web in blog posts articles, and such.

Submit an [issue] if you found a bug or have a great idea for a new feature!

## Development

Set up for local development:

1. Clone your `nf-ffperase` locally:

        git clone git@github.com:papaemmelab/nf-ffperase.git

1. Create a branch for local development:

        git pull
        git checkout -b name-of-your-bugfix-or-feature

    Now you can make your changes locally.

1. Create a test in:

        nf-ffperase/tests

1. Run [nf-test]:

        nf-test test

    For a dockerized version, verbose and with coverage metrics run:

        nf-test test --profile cloud --coverage --verbose:

1. Commit your changes and push your branch to GitHub (see our [`.gitmessage`] template):

        git add .
        git config commit.template .gitmessage
        git commit -m ":emoji-name: your short and nice description"
        git push origin name-of-your-bugfix-or-feature

    `emoji-name` should be one of the following:

    | emoji | name             | type of change              |
    | ----- | ---------------- | --------------------------- |
    | 🚀    | rocket           | new feature                 |
    | 🐛    | bug              | bug fix                     |
    | 📝    | memo             | changes to documentation    |
    | 🎨    | art              | formatting  no code change  |
    | 🔧    | wrench           | refactoring production code |
    | ✅    | white_check_mark | adding/editing test logic   |
    | 👕    | shirt            | no production code change   |
    | 💎    | gem              | bump to new version         |

    If you are suggesting a new version make sure you are following the [semantic versioning] guidelines and then update the [`VERSION`] file:

        git add papaemmelab/VERSION
        git commit -m ":gem: bump to version 0.1.0"

1. Submit a pull request through the GitHub website.

<!-- References -->
[`VERSION`]: ../VERSION
[`.gitmessage`]: ../.gitmessage
[semantic versioning]: http://semver.org/
[nf-test]: https://www.nf-test.com/
[issue]: https://github.com/papaemmelab/nf-ffperase/issues
