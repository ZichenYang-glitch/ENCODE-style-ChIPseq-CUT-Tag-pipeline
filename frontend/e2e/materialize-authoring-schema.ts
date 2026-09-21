// Cross-language smoke helper: use exactly the workbench's schema/default path.
import { readFileSync } from 'node:fs';

import {
  createDefaultObject,
  readWorkbenchSchema,
} from '../src/features/input-workbench/schemaContract';

const result = readWorkbenchSchema(JSON.parse(readFileSync(0, 'utf8')));
if (!result.ok) throw new Error(result.code);

process.stdout.write(JSON.stringify({
  config: createDefaultObject(result.value.configSchema),
  options: createDefaultObject(result.value.optionSchema),
}));
