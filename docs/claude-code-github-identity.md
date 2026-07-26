# A separate GitHub identity for Claude Code

**Question:** on a free GitHub account, can Claude Code act under its own profile instead of Marko's?

**Answer:** yes, two ways, both free. A *machine account* is explicitly carved out of the
one-free-account rule, and a *GitHub App* needs no second account at all.

## Why it currently shows as Marko

Two independent things decide attribution, which is why changing one is not enough:

| what | decided by | current value |
| --- | --- | --- |
| commit author | `git config user.name` / `user.email` | `Marko Tintor <marko@tintor.net>` |
| PR and comment author | the token `gh` is authenticated with | `tintor` |

So a commit can say one identity while the pull request that carries it says another. Both have to
be pointed at the new identity.

## Option 1: a machine account

GitHub's Terms of Service permit exactly this, as an exception to the one-account rule:

> One person or legal entity may maintain no more than one free Account (if you choose to control a
> machine account as well, that's fine, but it can only be used for running a machine).

and define it:

> A machine account is an Account set up by an individual human who accepts the Terms on behalf of
> the Account, provides a valid email address, and is responsible for its actions. A machine account
> is used exclusively for performing automated tasks.

Note the neighbouring sentence, which is why the account has to be registered by a human rather than
by the agent itself:

> Accounts registered by "bots" or other automated methods are not permitted. We do permit machine
> accounts.

Setup:

1. Register a second free account, e.g. `tintor-claude`, from a browser. It needs its own email
   address; a `+` alias such as `marko+claude@tintor.net` is usually accepted.
2. Add it to `tintor/algebra` as a collaborator with write access.
3. From that account, create a fine-grained personal access token scoped to `tintor/algebra` with
   Contents: read and write, Pull requests: read and write, Issues: read and write.
4. In the environment Claude Code runs in:

   ```sh
   gh auth login --with-token < token.txt      # or GH_TOKEN=... in the environment
   git config user.name  "tintor-claude"
   git config user.email "<the machine account's email>"
   ```

   Set `user.name`/`user.email` per repository rather than globally, so other work keeps Marko's
   identity.
5. Verify with `gh api user -q .login` and one throwaway commit and PR.

Trade-offs: commits, pull requests, comments and reviews all show the machine account, which is the
cleanest possible separation. But the account counts as a collaborator, its token has to be stored
somewhere, and it cannot approve pull requests in a way that satisfies a required-review rule if the
rule excludes the author — the same as any other collaborator.

## Option 2: a GitHub App

An App is not an account, so there is nothing to register under the ToS account rules, and Apps can
be created from a personal account. Authenticating as an *installation* attributes the work to the
App:

> API requests made by an app installation are attributed to the app.

> the timeline of the issue would state that your app updated the status

This is the mechanism Dependabot and similar tools use, and the actor appears in the UI as the App
rather than as a person.

Setup sketch: create the App under Settings → Developer settings → GitHub Apps, give it Contents,
Pull requests and Issues write permissions, install it on `tintor/algebra`, then have the runtime
exchange the App's private key for a short-lived installation token and use that as `GH_TOKEN`.
Commits still need `user.name`/`user.email` set separately, since commit authorship is not decided by
the token.

Trade-offs: no second account and no long-lived secret, tokens expire after an hour, and permissions
are scoped per repository. The cost is the token-exchange step, which needs a small script and the
App's private key.

## Recommendation

For this repository, **option 1** is the lower-effort choice: one browser registration, one token,
two `git config` lines, and every trace of automated work then reads `tintor-claude`. Option 2 is
better if the credential handling matters more than the setup time, since an installation token is
short lived and cannot be reused elsewhere.

Either way, keep `user.name`/`user.email` set locally in this repository, so the identity change does
not leak into unrelated work.

## Verifying it worked

```sh
gh api user -q .login          # the new identity
git log -1 --format='%an <%ae>'
gh pr view <n> --json author -q .author.login
```

All three have to show the new identity; if only some do, one of the two mechanisms above was
missed.
