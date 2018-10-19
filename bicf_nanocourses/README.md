## BICF Nanocourse Pages

Learn more about GitLab Pages at https://pages.gitlab.io and the official
documentation https://docs.gitlab.com/ce/user/project/pages/.


## Structure

`mkdocs.yml` is a configuration file for the mkdocs static site generator. In here
you will need to set your site's title, and the URL for the site (which you will set) in
the GitLab project settings.

`docs` is a folder containing your markdown pages that you wish to publish to the web. The
front page for your site should be name `index.md`.


## GitLab CI

This project's static Pages are built by [GitLab CI][ci], following the steps
defined in [`.gitlab-ci.yml`](.gitlab-ci.yml):


## Building locally

To work locally with this project on a BioHPC machine, you'll have to follow the steps below:

1. Fork, clone or download this project
1. Load the anaconda python module `module add python/2.7.x-anaconda`
1. Install the python mkdocs package to your home dir `pip install --user mkdocs`
1. Preview your project: `mkdocs serve`,
   your site can be accessed under `localhost:8000`
1. Add content
1. Generate the website: `mkdocs build` (optional)

Read more at MkDocs [documentation][].


## Did you fork this project?

If you forked this project for your own use, please go to your project's
**Settings** and remove the forking relationship, which won't be necessary
unless you want to contribute back to the upstream project.

## Troubleshooting

1. CSS is missing! That means two things:

    Either that you have wrongly set up the CSS URL in your templates, or
    your static generator has a configuration option that needs to be explicitly
    set in order to serve static assets under a relative URL.

[ci]: https://about.gitlab.com/gitlab-ci/
[mkdocs]: http://www.mkdocs.org
[install]: http://www.mkdocs.org/#installation
[documentation]: http://www.mkdocs.org
[userpages]: https://docs.gitlab.com/ce/user/project/pages/introduction.html#user-or-group-pages
[projpages]: https://docs.gitlab.com/ce/user/project/pages/introduction.html#project-pages

---

Forked from https://gitlab.com/morph027/mkdocs
