# DoppelBot: Fine-tune an LLM to replace your CEO

*(quick links:
[add to your own Slack](https://github.com/modal-labs/doppel-bot#usage);
[source code](https://github.com/modal-labs/doppel-bot))*

Internally at Modal, we spend a *lot* of time talking to each other on Slack.
Now, with the advent of open-source large language models, we had started to
wonder if all of this wasn't a bit redundant. Could we have these language
models bike-shed on Slack for us, so we could spend our time on higher leverage
activities such as
[paddleboarding in Tahiti](https://x.com/modal/status/1642262543757352960)
instead?

To test this, we fine-tuned
[Llama 3.1](https://ai.meta.com/blog/meta-llama-3-1/) on
[Erik](https://twitter.com/bernhardsson)'s Slack messages, and `@erik-bot` was
born.

![erik-bot](https://modal-cdn.com/erik-bot-1.jpeg)

Since then, `@erik-bot` has been an invaluable asset to us, in areas ranging
from [API design](https://modal-cdn.com/erik-bot-2.png) to
[legal advice](https://modal-cdn.com/erik-bot-3.png) to thought leadership.

![erik-bot-3](https://modal-cdn.com/erik-bot-4.png)

We were planning on releasing the weights for `@erik-bot` to the world, but all
our metrics have been going up and to the right a little too much since we've
launched him...

So, we are releasing the next best thing. `DoppelBot` is a Slack bot that you
can install in your own workspace, and fine-tune on your own Slack messages.
Follow the instructions [here](https://github.com/modal-labs/doppel-bot#usage)
to replace your own CEO with an LLM today.

All the components—scraping, fine-tuning, inference and slack event handlers run
on Modal, and the code itself is open-source and available
[here](https://github.com/modal-labs/doppel-bot). If you're new to Modal, it's
worth reiterating that **all of these components are also serverless and scale
to zero**. This means that you can deploy and forget about them, because you'll
only pay for compute when your app is used!

## How it works

DoppelBot uses the Slack SDK to scrape messages from a Slack workspace, and
converts them into prompt/response pairs. It uses these to fine-tune a language
model using [Low-Rank Adaptation (LoRA)](https://arxiv.org/abs/2106.09685), a
technique that produces a small adapter that can be merged with the base model
when needed, instead of modifying all the parameters in the base model. The
fine-tuned adapters for each user are stored in a Modal
[Volume](/docs/guide/volumes). When a user `@`s the bot,
Slack sends a webhook call to Modal, which loads the adapter for that user and
generates a response.

We go into detail into each of these steps below, and provide commands for
running each of them individually. To follow along,
[clone the repo](https://github.com/modal-labs/doppel-bot) and
[set up a Slack token](https://github.com/modal-labs/doppel-bot#create-a-slack-app)
for yourself.

### Scraping slack

<GuideGithubLink url="https://github.com/modal-labs/doppel-bot/blob/main/src/scrape.py" />

The scraper uses Modal's [`.map()`](/docs/guide/scale#scaling-out) to fetch
messages from all public channels in parallel. Each thread is split into
contiguous messages from the target users and continguous messages from other
users. These will be fed into the model as prompts in the following format:

```
[system]: You are {user}, employee at a fast-growing startup. Below is an input conversation that takes place in the company's internal Slack. Write a response that appropriately continues the conversation.

[user]: <slack thread>

[assistant]: <target user's response>
```

Initial versions of the model were prone to generating short responses
— unsurprising, because a majority of Slack communication is pretty terse.
Adding a minimum character length for the target user's messages fixed this.

If you're following along at home, you can run the scraper with the following
command:

```bash
modal run -m src.scrape::scrape --user="<user>"
```

Scraped results are stored in a Modal
[Volume](/docs/guide/volumes), so they can be used by the next step.

### Fine-tuning

<GuideGithubLink url="https://github.com/modal-labs/doppel-bot/blob/main/src/finetune.py" />

Next, we use the prompts to fine-tune a language model. We chose
[Llama 3.1](https://ai.meta.com/blog/meta-llama-3-1/) because of its permissive license and high quality relative to its small size. Fine-tuning is
done using [Low-Rank Adaptation (LoRA)](https://arxiv.org/abs/2106.09685), a
[parameter-efficient fine-tuning](https://huggingface.co/blog/peft) technique
that produces a small adapter that can be merged with the base model when needed
(~60MB for the rank we're using).

Our fine-tuning implementation uses [torchtune](https://github.com/pytorch/torchtune), a new PyTorch library for easily configuring fine-tuning runs.

Because of the typically small sample sizes we're working with, training for
longer than a couple hundred steps (with our batch size of 128) quickly led to
overfitting. Admittedly, we haven't thoroughly evaluated the hyperparameter
space yet — do reach out to us if you're interested in collaborating on this!

![train-loss](../../assets/docs/train-loss.png)

To try this step yourself, run:

```bash
modal run -m src.finetune --user="<user>"
```

### Inference

<GuideGithubLink url="https://github.com/modal-labs/doppel-bot/blob/main/src/inference.py" />

We use [vLLM](https://github.com/vllm-project/vllm) as our inference engine, which now comes with support for dynamically swapping LoRA adapters [out of the box](https://docs.vllm.ai/en/latest/features/lora.html).

With parametrized functions, every user model gets its own pool of containers
that scales up when there are incoming requests, and scales to 0 when there's
none. Here's what that looks like stripped down to the essentials:

```python notest
@app.cls(gpu="L40S")
class Model():
    @modal.enter()
    def enter(self):
        self.engine = AsyncLLMEngine.from_engine_args(AsyncEngineArgs(...))
        self.loras: dict[str, int] = dict()  # per replica LoRA identifier

    @modal.method()
    def generate(self, input: str):
        if (ident := f"{user}-{team_id}") not in self.loras:
            self.loras[ident] = len(self.loras) + 1

        lora_request = LoRARequest(
            ident, self.loras[ident], lora_local_path=checkpoint_path
        )

        tokenizer = await self.engine.get_tokenizer(lora_request=lora_request)

        prompt = tokenizer.apply_chat_template(
            conversation=inpt, tokenize=False, add_generation_prompt=True
        )

        results_generator = self.engine.generate(prompt, lora_request=lora_request,)
```

If you've fine-tuned a model already in the previous step, you can run inference
using it now:

```bash
modal run -m src.inference --user="<user>"
```

(We have a list of sample inputs in the file, but you can also try it out with
your own messages!)

### Slack Bot

<GuideGithubLink url="https://github.com/modal-labs/doppel-bot/blob/main/src/bot.py" />

Finally, it all comes together in
[`bot.py`](https://github.com/modal-labs/doppel-bot/blob/main/src/bot.py). As
you might have guessed, all events from Slack are handled by serverless Modal
functions. We handle 3 types of events:

* [`url_verification`](https://github.com/modal-labs/doppel-bot/blob/24609583c43c0e722f56f85a1c00bb55b46c7754/src/bot.py#L112):
  To verify that this is a Slack app, Slack expects us to return a challenge
  string.
* [`app_mention`](https://github.com/modal-labs/doppel-bot/blob/main/src/bot.py#L118):
  When the bot is mentioned in a channel, we retrieve the recent messages from
  that thread, do some basic cleaning and call the user's model to generate a
  response.

```python notest
model = OpenLlamaModel.remote(user, team_id)
result = model.generate(messages)
```

* [`doppel` slash command](https://github.com/modal-labs/doppel-bot/blob/main/src/bot.py#L182):
  This command kicks off the scraping -> finetuning pipeline for the user.

To deploy the slackbot in its entirety, you need to run:

```shell
modal deploy -m src.bot
```

<div>

### Multi-Workspace Support

</div>

Everything we've talked about so far is for a single-workspace Slack app. To
make it work with multiple workspaces, we'll need to handle
[workspace installation and authentication with OAuth](https://api.slack.com/authentication/oauth-v2),
and also store some state for each workspace.

Luckily, Slack's [Bolt](https://slack.dev/bolt-python/concepts) framework
provides a complete (but frugally documented) OAuth implemention. A neat feature
is that the OAuth state can be backed by a file system, so all we need to do is
[point Bolt](https://github.com/modal-labs/doppel-bot/blob/24609583c43c0e722f56f85a1c00bb55b46c7754/src/bot.py#L78)
at a Modal [Volume](/docs/guide/volumes), and then we don't need to worry about
managing this state ourselves.

To store state for each workspace, we're using [Neon](https://neon.tech/), a
serverless Postgres database that's really easy to set up and *just works*. If
you're interested in developing a multi-workspace app,
[follow our instructions](https://github.com/modal-labs/doppel-bot#optional-multi-workspace-app)
on how to set up Neon with Modal.

## Next Steps

If you've made it this far, you have just found a way to increase your team's
productivity by 10x! Congratulations on the well-earned vacation! 🎉

If you're interested in learning more about Modal, check out our [docs](/docs)
and other [examples](/examples).
